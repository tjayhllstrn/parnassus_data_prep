#ifndef CutStudy_h
#define CutStudy_h

// ---------------------------------------------------------------------------
// CutStudy -- machinery for measuring what each particle-level cut does.
//
// CutManager::filter_particles() only hands back the survivors, so there is no
// way to ask "which cut emptied this event?". This header provides the plumbing
// to evaluate cuts one at a time and tally the result; the list of cuts to
// study lives in the driver macro (macros/cut_study.C), not here.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "TSystem.h"

#include "CutManager.h"
#include "Structs.h"

// ---------------------------------------------------------------------------
// CutManager deliberately keeps its individual cuts `protected` so production
// code can only go through filter_particles(). The study needs to call them
// one at a time, so re-expose them here. This changes nothing about how the
// cuts behave -- the run period / torus bending state is still set up the
// normal way, with set_run().
// ---------------------------------------------------------------------------
class ProbeableCutManager : public CutManager {
 public:
  using CutManager::is_in_FD;

  // cuts applied to every particle regardless of PID
  using CutManager::FTLeptonPhotonCut;
  using CutManager::maxPCut;
  using CutManager::VrCut_outliers;

  // whole per-PID cut groups
  using CutManager::apply_electron_cuts;
  using CutManager::apply_pion_cuts;
  using CutManager::apply_charged_hadron_cuts;
  using CutManager::apply_photon_cuts;

  // the individual cuts those groups are built from
  using CutManager::DC_fiducial_cut;
  using CutManager::EleSampFrac_pass2;
  using CutManager::minEpcal;
  using CutManager::caloEdges;
  using CutManager::chi2pid;
  using CutManager::VzCut_outliers;
  using CutManager::betaCut_outlier;

  // defined but NOT currently applied by production -- available for what-if
  using CutManager::thetaRangeCut;
  using CutManager::VzCut;
  using CutManager::vz_e_pi;
  using CutManager::photonBetaCut;
  using CutManager::leptonBetaCut;
  using CutManager::photonElectronAngle;
  using CutManager::minNeutronEnergy;
  using CutManager::minEleMomentum;
  using CutManager::minPiMomentum;
  using CutManager::Ele3calo;
  using CutManager::photonMinEtot;
  using CutManager::hadronStatus;
};

// ---------------------------------------------------------------------------
// One cut under study.
//   pass    : returns true when the particle is KEPT
//   pids    : PIDs this cut is meaningful for; empty means "every particle"
//   fd_only : only evaluate for particles in the forward detector, mirroring
//             the `if (is_in_FD(particle))` guards inside CutManager
//
//   in_stack: part of the production cut stack. The in_stack probes together
//             must reproduce filter_particles(); process_event() checks this
//             every event. Composite and what-if probes are not in the stack.
//
// A particle a probe does not apply to is counted as kept -- which is what
// filter_particles() effectively does. This gating matters: several CutManager
// helpers return false for particles they were never meant to judge (chi2pid()
// returns false for a non-pion without a chi2 value, DC_fiducial_cut_edge()
// returns false for an unlisted PID), so calling them unguarded would report
// wildly inflated rejection rates.
// ---------------------------------------------------------------------------
struct CutProbe {
  // (cut manager, particle, scattered electron, was an electron found?)
  using Fn = std::function<bool(ProbeableCutManager&, const part&, const part&, bool)>;

  std::string      name;
  std::string      group;
  std::vector<int> pids;
  bool             fd_only = false;
  Fn               pass;
  bool             in_stack = true;

  bool applies(ProbeableCutManager& cm, const part& particle) const {
    if (fd_only && !cm.is_in_FD(particle)) return false;
    if (pids.empty()) return true;
    return std::find(pids.begin(), pids.end(), particle.pid) != pids.end();
  }
};

// ---------------------------------------------------------------------------
// Accumulates, over many events, what each probe does
//   ALONE : applied to the uncut particles with every other cut off
//   N-1   : applied to the particles that pass every OTHER production cut, i.e.
//           what the final dataset loses because of this cut (for a probe not
//           in the stack: what it would lose if the cut were added)
// ---------------------------------------------------------------------------
class CutStudy {
 public:
  explicit CutStudy(std::vector<CutProbe> probes) : _probes(std::move(probes)) {
    _stats.resize(_probes.size());
  }

  // Choose the combination of cuts whose emptied events process_event() should
  // flag, so the caller can save them. `spec` is a comma-separated list of
  // selectors, each matching either a whole group or one cut as "group/name":
  //   "ALL"                  -> the full production stack (filter_particles)
  //   "global,vertex"        -> every probe in those two groups, ANDed together
  //   "photon/caloEdges"     -> that one cut on its own
  //   ""                     -> capture disabled
  // Returns false if a selector matched no probe (nothing is captured then).
  bool set_capture_cuts(const std::string& spec) {
    _capture_full_stack = false;
    _capture_probes.clear();
    bool ok = true;
    std::stringstream ss(spec);
    std::string token;
    while (std::getline(ss, token, ',')) {
      // trim surrounding whitespace
      const size_t b = token.find_first_not_of(" \t");
      if (b == std::string::npos) continue;
      token = token.substr(b, token.find_last_not_of(" \t") - b + 1);

      if (token == "ALL") { _capture_full_stack = true; continue; }
      size_t n_matched = 0;
      for (size_t ic = 0; ic < _probes.size(); ic++) {
        if (_probes[ic].group == token ||
            _probes[ic].group + "/" + _probes[ic].name == token) {
          _capture_probes.push_back(ic);
          n_matched++;
        }
      }
      if (n_matched == 0) {
        std::cerr << "  [CutStudy] capture selector '" << token
                  << "' matched no cut in the registry" << std::endl;
        ok = false;
      }
    }
    return ok;
  }

  bool capture_enabled() const { return _capture_full_stack || !_capture_probes.empty(); }

  // Feed one event's *unfiltered* reco particles, with is_scattered_electron
  // already tagged the same way macros/hipo2tree_parnassus.C tags it.
  // Returns true when the capture combination (set_capture_cuts) leaves this
  // event with no surviving particle -- i.e. the event is one to save.
  bool process_event(ProbeableCutManager& cm, const std::vector<part>& particles) {
    _n_events++;
    _n_particles += (long)particles.size();
    if (particles.empty()) _n_events_empty_nocuts++;

    // Reference point: the real production cut stack, unmodified.
    const long n_production = (long)cm.filter_particles(particles).size();
    const bool empty_baseline = n_production == 0;
    if (empty_baseline) _n_events_empty_baseline++;
    _n_kept_production += n_production;

    // Reproduce how filter_particles decides whether the electron-dependent
    // cuts can be applied at all (CutManager.C, top of filter_particles).
    bool electron_found = !(particles.empty() || particles[0].is_scattered_electron == -1);
    part electron;
    if (electron_found) {
      for (const auto& particle : particles)
        if (particle.is_scattered_electron == 1) electron = particle;
    }

    // Evaluate every probe once per particle; the ALONE, N-1 and stack numbers
    // all come from the same verdicts.
    const size_t np = _probes.size();
    std::vector<long> n_survivors(np, 0);
    std::vector<char> judged(np), kept(np);
    long n_kept_stack = 0;
    for (const auto& particle : particles) {
      int n_stack_fail = 0;
      for (size_t ic = 0; ic < np; ic++) {
        const CutProbe& probe = _probes[ic];
        judged[ic] = probe.applies(cm, particle);
        kept[ic]   = !judged[ic] || probe.pass(cm, particle, electron, electron_found); // not judged -> kept
        if (kept[ic]) n_survivors[ic]++;
        if (judged[ic]) {
          ProbeStats& st = _stats[ic];
          st.n_seen++;
          st.seen_by_pid[particle.pid]++;
          if (!kept[ic]) { st.n_removed++; st.removed_by_pid[particle.pid]++; }
        }
        if (probe.in_stack && !kept[ic]) n_stack_fail++;
      }
      if (n_stack_fail == 0) n_kept_stack++;

      // N-1: only particles that pass every production cut other than this one
      for (size_t ic = 0; ic < np; ic++) {
        if (!judged[ic]) continue;
        const int others_fail = n_stack_fail - ((_probes[ic].in_stack && !kept[ic]) ? 1 : 0);
        if (others_fail != 0) continue;
        ProbeStats& st = _stats[ic];
        st.n1_base++;
        st.n1_base_by_pid[particle.pid]++;
        if (!kept[ic]) { st.n1_removed++; st.n1_removed_by_pid[particle.pid]++; }
      }
    }
    for (size_t ic = 0; ic < np; ic++)
      if (n_survivors[ic] == 0) _stats[ic].n_events_empty++;
    _n_kept_stack += n_kept_stack;
    if (n_kept_stack != n_production) _n_events_stack_mismatch++;

    // Does the chosen capture combination empty this event? For "ALL" the
    // answer is the baseline we already computed, so it costs nothing.
    bool capture = false;
    if (_capture_full_stack) {
      capture = empty_baseline;
    } else if (!_capture_probes.empty()) {
      capture = true;
      for (const auto& particle : particles) {
        bool keep = true;
        for (size_t ic : _capture_probes) {
          const CutProbe& probe = _probes[ic];
          if (!probe.applies(cm, particle)) continue;   // not judged -> kept
          if (!probe.pass(cm, particle, electron, electron_found)) { keep = false; break; }
        }
        if (keep) { capture = false; break; }  // one survivor is enough, stop early
      }
    }
    if (capture) _n_events_captured++;
    return capture;
  }

  void report(std::ostream& os) const {
    auto pct = [](long a, long b){ return b > 0 ? 100.0 * a / b : 0.0; };

    const int W = 136;
    os << "\n" << std::string(W, '=') << "\n";
    os << "  Particle-level cut study -- each cut ALONE (every other cut off) and N-1 (every other cut on)\n";
    os << std::string(W, '=') << "\n";
    os << "  events processed                  : " << _n_events << "\n";
    os << "  reco particles loaded (pre-cut)   : " << _n_particles << "\n";
    os << "  reco particles kept, ALL cuts     : " << _n_kept_production
       << "  (" << std::fixed << std::setprecision(2) << pct(_n_kept_production, _n_particles)
       << "%, from filter_particles)\n";
    os << "  registry reproduces ALL cuts      : "
       << (_n_events_stack_mismatch == 0 ? "yes" : "NO")
       << "  (registry stack kept " << _n_kept_stack << "; events that differ: "
       << _n_events_stack_mismatch << ")\n";
    if (_n_events_stack_mismatch != 0)
      os << "  WARNING: the registry no longer matches filter_particles; the N-1 numbers are not trustworthy\n";
    os << "  events with 0 particles, NO cuts  : " << _n_events_empty_nocuts
       << "  (" << pct(_n_events_empty_nocuts, _n_events) << "%)\n";
    os << "  events with npflow==0, ALL cuts   : " << _n_events_empty_baseline
       << "  (" << pct(_n_events_empty_baseline, _n_events) << "%)   <- production baseline\n";
    if (capture_enabled()) {
      os << "  events captured to the ROOT file  : " << _n_events_captured
         << "  (" << pct(_n_events_captured, _n_events) << "%)\n";
    }
    os << std::string(W, '-') << "\n";
    os << "  ALONE: 'npflow==0' counts events this cut alone leaves with no particle; 'seen' = particles the\n";
    os << "         cut judged (after PID / forward-detector gating), 'removed' = of those, rejected.\n";
    os << "  N-1:   'base' = judged particles that pass every OTHER production cut, 'removed' = of those,\n";
    os << "         rejected by this cut, '% final' = removed / particles kept by ALL cuts. For a production\n";
    os << "         cut this is what the final dataset loses to it; for a composite or what-if probe\n";
    os << "         (marked *) it is what adding it on top of the production cuts would remove.\n";
    os << std::string(W, '-') << "\n";
    os << "  " << std::left << std::setw(10) << "group"
       << std::setw(28) << "cut"
       << std::right
       << std::setw(11) << "npflow==0"
       << std::setw(8)  << "% ev"
       << std::setw(12) << "seen"
       << std::setw(11) << "removed"
       << std::setw(8)  << "% rm"
       << "  |"
       << std::setw(11) << "N-1 base"
       << std::setw(11) << "removed"
       << std::setw(8)  << "% rm"
       << std::setw(9)  << "% final" << "\n";

    std::string last_group;
    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const CutProbe&   probe = _probes[ic];
      const ProbeStats& st    = _stats[ic];
      if (probe.group != last_group) {
        os << "  " << std::string(W - 4, '.') << "\n";
        last_group = probe.group;
      }
      os << "  " << std::left << std::setw(10) << probe.group
         << std::setw(28) << (probe.in_stack ? probe.name : probe.name + " *")
         << std::right
         << std::setw(11) << st.n_events_empty
         << std::setw(8)  << pct(st.n_events_empty, _n_events)
         << std::setw(12) << st.n_seen
         << std::setw(11) << st.n_removed
         << std::setw(8)  << pct(st.n_removed, st.n_seen)
         << "  |"
         << std::setw(11) << st.n1_base
         << std::setw(11) << st.n1_removed
         << std::setw(8)  << pct(st.n1_removed, st.n1_base)
         << std::setw(9)  << pct(st.n1_removed, _n_kept_production) << "\n";
    }

    os << std::string(W, '-') << "\n";
    os << "  Particles removed per PID, ALONE and N-1 (cuts that removed nothing are omitted):\n";
    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const CutProbe&   probe = _probes[ic];
      const ProbeStats& st    = _stats[ic];
      if (st.n_removed == 0) continue;
      os << "    " << probe.group << " / " << probe.name << (probe.in_stack ? "" : " *") << "\n";
      for (const auto& kv : st.seen_by_pid) {
        const long removed    = lookup(st.removed_by_pid, kv.first);
        const long n1_base    = lookup(st.n1_base_by_pid, kv.first);
        const long n1_removed = lookup(st.n1_removed_by_pid, kv.first);
        if (removed == 0 && n1_removed == 0) continue;
        os << "      pid " << std::setw(6) << kv.first
           << " : alone " << std::setw(10) << removed << " of " << std::setw(10) << kv.second
           << " (" << std::setw(6) << pct(removed, kv.second) << "%)"
           << "   N-1 " << std::setw(10) << n1_removed << " of " << std::setw(10) << n1_base
           << " (" << std::setw(6) << pct(n1_removed, n1_base) << "%)\n";
      }
    }
    os << std::string(W, '=') << "\n";
  }

  // summary.csv    : one row per cut (plus two reference rows)
  // removed_by_pid.csv : one row per (cut, pid)
  void write_csv(const std::string& dir) const {
    gSystem->mkdir(dir.c_str(), true);

    const std::string summary_path = dir + "/summary.csv";
    std::ofstream summary(summary_path.c_str());
    summary << "group,cut,in_stack,n_events,n_events_npflow0,frac_npflow0,"
            << "n_particles_seen,n_particles_removed,frac_removed,"
            << "n1_base,n1_removed,frac_n1_removed,frac_of_final\n";
    summary << "reference,no_cuts,0," << _n_events << "," << _n_events_empty_nocuts << ","
            << frac(_n_events_empty_nocuts, _n_events) << "," << _n_particles << ",0,0,0,0,0,0\n";
    summary << "reference,all_cuts_baseline,0," << _n_events << "," << _n_events_empty_baseline << ","
            << frac(_n_events_empty_baseline, _n_events) << "," << _n_particles << ","
            << (_n_particles - _n_kept_production) << ","
            << frac(_n_particles - _n_kept_production, _n_particles) << ",0,0,0,0\n";
    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const ProbeStats& st = _stats[ic];
      summary << _probes[ic].group << "," << _probes[ic].name << "," << (_probes[ic].in_stack ? 1 : 0) << ","
              << _n_events << "," << st.n_events_empty << ","
              << frac(st.n_events_empty, _n_events) << ","
              << st.n_seen << "," << st.n_removed << ","
              << frac(st.n_removed, st.n_seen) << ","
              << st.n1_base << "," << st.n1_removed << ","
              << frac(st.n1_removed, st.n1_base) << ","
              << frac(st.n1_removed, _n_kept_production) << "\n";
    }
    summary.close();

    const std::string pid_path = dir + "/removed_by_pid.csv";
    std::ofstream bypid(pid_path.c_str());
    bypid << "group,cut,pid,n_seen,n_removed,frac_removed,n1_base,n1_removed,frac_n1_removed\n";
    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const ProbeStats& st = _stats[ic];
      for (const auto& kv : st.seen_by_pid) {
        const long removed    = lookup(st.removed_by_pid, kv.first);
        const long n1_base    = lookup(st.n1_base_by_pid, kv.first);
        const long n1_removed = lookup(st.n1_removed_by_pid, kv.first);
        bypid << _probes[ic].group << "," << _probes[ic].name << "," << kv.first << ","
              << kv.second << "," << removed << "," << frac(removed, kv.second) << ","
              << n1_base << "," << n1_removed << "," << frac(n1_removed, n1_base) << "\n";
      }
    }
    bypid.close();

    std::cout << "  wrote " << summary_path << "\n";
    std::cout << "  wrote " << pid_path << "\n";
  }

  long n_events()                const { return _n_events; }
  long n_events_empty_baseline() const { return _n_events_empty_baseline; }
  long n_events_captured()       const { return _n_events_captured; }

 private:
  struct ProbeStats {
    long n_events_empty = 0;  // events this cut alone left with zero particles
    long n_seen         = 0;  // particles the cut actually judged
    long n_removed      = 0;  // of those, how many it rejected
    long n1_base        = 0;  // judged particles that pass every other production cut
    long n1_removed     = 0;  // of those, how many this cut rejected
    std::map<int,long> seen_by_pid;
    std::map<int,long> removed_by_pid;
    std::map<int,long> n1_base_by_pid;
    std::map<int,long> n1_removed_by_pid;
  };

  static double frac(long a, long b) { return b > 0 ? (double)a / (double)b : 0.0; }
  static long lookup(const std::map<int,long>& m, int key) {
    auto it = m.find(key);
    return it == m.end() ? 0 : it->second;
  }

  std::vector<CutProbe>   _probes;
  std::vector<ProbeStats> _stats;
  long _n_events                = 0;
  long _n_particles             = 0;
  long _n_events_empty_nocuts   = 0;
  long _n_events_empty_baseline = 0;
  long _n_kept_production       = 0;  // particles kept by filter_particles
  long _n_kept_stack            = 0;  // particles kept by the registry's in_stack probes
  long _n_events_stack_mismatch = 0;  // events where those two differ

  // capture selection (see set_capture_cuts)
  bool                _capture_full_stack = false;
  std::vector<size_t> _capture_probes;
  long                _n_events_captured  = 0;
};

#endif
