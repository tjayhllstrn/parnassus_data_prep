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
  using CutManager::thetaRangeCut;
  using CutManager::FTLeptonPhotonCut;
  using CutManager::maxPtCut;

  // whole per-PID cut groups
  using CutManager::apply_electron_cuts;
  using CutManager::apply_pion_cuts;
  using CutManager::apply_photon_cuts;
  using CutManager::apply_neutron_cuts;

  // the individual cuts those groups are built from
  using CutManager::DC_fiducial_cut;
  using CutManager::EleSampFrac_pass2;
  using CutManager::VzCut;
  using CutManager::minEpcal;
  using CutManager::caloEdges;
  using CutManager::chi2pid;
  using CutManager::vz_e_pi;
  using CutManager::photonBetaCut;
  using CutManager::photonElectronAngle;
  using CutManager::minNeutronEnergy;

  // defined but NOT currently applied by production -- available for what-if
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
// A particle a probe does not apply to is counted as kept -- which is what
// filter_particles() effectively does. This gating matters: several CutManager
// helpers return false for particles they were never meant to judge (chi2pid()
// returns false for anything that is not a charged pion, DC_fiducial_cut_XY()
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

  bool applies(ProbeableCutManager& cm, const part& particle) const {
    if (fd_only && !cm.is_in_FD(particle)) return false;
    if (pids.empty()) return true;
    return std::find(pids.begin(), pids.end(), particle.pid) != pids.end();
  }
};

// ---------------------------------------------------------------------------
// Accumulates, over many events, what each probe does when applied ALONE.
// ---------------------------------------------------------------------------
class CutStudy {
 public:
  explicit CutStudy(std::vector<CutProbe> probes) : _probes(std::move(probes)) {
    _stats.resize(_probes.size());
  }

  // Feed one event's *unfiltered* reco particles, with is_scattered_electron
  // already tagged the same way macros/hipo2tree_parnassus.C tags it.
  void process_event(ProbeableCutManager& cm, const std::vector<part>& particles) {
    _n_events++;
    _n_particles += (long)particles.size();
    if (particles.empty()) _n_events_empty_nocuts++;

    // Reference point: the real production cut stack, unmodified.
    if (cm.filter_particles(particles).empty()) _n_events_empty_baseline++;

    // Reproduce how filter_particles decides whether the electron-dependent
    // cuts can be applied at all (CutManager.C, top of filter_particles).
    bool electron_found = !(particles.empty() || particles[0].is_scattered_electron == -1);
    part electron;
    if (electron_found) {
      for (const auto& particle : particles)
        if (particle.is_scattered_electron == 1) electron = particle;
    }

    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const CutProbe& probe = _probes[ic];
      ProbeStats&     st    = _stats[ic];

      long n_survivors = 0;
      for (const auto& particle : particles) {
        if (!probe.applies(cm, particle)) { n_survivors++; continue; } // not judged -> kept
        st.n_seen++;
        st.seen_by_pid[particle.pid]++;
        if (probe.pass(cm, particle, electron, electron_found)) {
          n_survivors++;
        } else {
          st.n_removed++;
          st.removed_by_pid[particle.pid]++;
        }
      }
      if (n_survivors == 0) st.n_events_empty++;
    }
  }

  void report(std::ostream& os) const {
    auto pct = [](long a, long b){ return b > 0 ? 100.0 * a / b : 0.0; };

    os << "\n" << std::string(94, '=') << "\n";
    os << "  Particle-level cut study -- each cut applied ALONE (every other cut off)\n";
    os << std::string(94, '=') << "\n";
    os << "  events processed                  : " << _n_events << "\n";
    os << "  reco particles loaded (pre-cut)   : " << _n_particles << "\n";
    os << "  events with 0 particles, NO cuts  : " << _n_events_empty_nocuts
       << "  (" << std::fixed << std::setprecision(2)
       << pct(_n_events_empty_nocuts, _n_events) << "%)\n";
    os << "  events with npflow==0, ALL cuts   : " << _n_events_empty_baseline
       << "  (" << pct(_n_events_empty_baseline, _n_events) << "%)   <- production baseline\n";
    os << std::string(94, '-') << "\n";
    os << "  Per-cut, applied alone. 'npflow==0' counts events left with no surviving particle;\n";
    os << "  it can never be below the NO-cuts number or above the ALL-cuts baseline.\n";
    os << "  'seen' = particles the cut actually judged (after PID / forward-detector gating).\n";
    os << std::string(94, '-') << "\n";
    os << "  " << std::left << std::setw(10) << "group"
       << std::setw(26) << "cut"
       << std::right
       << std::setw(11) << "npflow==0"
       << std::setw(9)  << "% ev"
       << std::setw(12) << "seen"
       << std::setw(11) << "removed"
       << std::setw(9)  << "% rm" << "\n";

    std::string last_group;
    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const CutProbe&   probe = _probes[ic];
      const ProbeStats& st    = _stats[ic];
      if (probe.group != last_group) {
        os << "  " << std::string(90, '.') << "\n";
        last_group = probe.group;
      }
      os << "  " << std::left << std::setw(10) << probe.group
         << std::setw(26) << probe.name
         << std::right
         << std::setw(11) << st.n_events_empty
         << std::setw(9)  << pct(st.n_events_empty, _n_events)
         << std::setw(12) << st.n_seen
         << std::setw(11) << st.n_removed
         << std::setw(9)  << pct(st.n_removed, st.n_seen) << "\n";
    }

    os << std::string(94, '-') << "\n";
    os << "  Particles removed per PID (cuts that removed nothing are omitted):\n";
    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const CutProbe&   probe = _probes[ic];
      const ProbeStats& st    = _stats[ic];
      if (st.n_removed == 0) continue;
      os << "    " << probe.group << " / " << probe.name << "\n";
      for (const auto& kv : st.removed_by_pid) {
        long seen = 0;
        auto it = st.seen_by_pid.find(kv.first);
        if (it != st.seen_by_pid.end()) seen = it->second;
        os << "      pid " << std::setw(6) << kv.first
           << " : removed " << std::setw(10) << kv.second
           << " of " << std::setw(10) << seen
           << "  (" << pct(kv.second, seen) << "%)\n";
      }
    }
    os << std::string(94, '=') << "\n";
  }

  // summary.csv    : one row per cut (plus two reference rows)
  // removed_by_pid.csv : one row per (cut, pid)
  void write_csv(const std::string& dir) const {
    gSystem->mkdir(dir.c_str(), true);

    const std::string summary_path = dir + "/summary.csv";
    std::ofstream summary(summary_path.c_str());
    summary << "group,cut,n_events,n_events_npflow0,frac_npflow0,"
            << "n_particles_seen,n_particles_removed,frac_removed\n";
    summary << "reference,no_cuts," << _n_events << "," << _n_events_empty_nocuts << ","
            << frac(_n_events_empty_nocuts, _n_events) << ",0,0,0\n";
    summary << "reference,all_cuts_baseline," << _n_events << "," << _n_events_empty_baseline << ","
            << frac(_n_events_empty_baseline, _n_events) << ",0,0,0\n";
    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const ProbeStats& st = _stats[ic];
      summary << _probes[ic].group << "," << _probes[ic].name << ","
              << _n_events << "," << st.n_events_empty << ","
              << frac(st.n_events_empty, _n_events) << ","
              << st.n_seen << "," << st.n_removed << ","
              << frac(st.n_removed, st.n_seen) << "\n";
    }
    summary.close();

    const std::string pid_path = dir + "/removed_by_pid.csv";
    std::ofstream bypid(pid_path.c_str());
    bypid << "group,cut,pid,n_seen,n_removed,frac_removed\n";
    for (size_t ic = 0; ic < _probes.size(); ic++) {
      const ProbeStats& st = _stats[ic];
      for (const auto& kv : st.seen_by_pid) {
        long removed = 0;
        auto it = st.removed_by_pid.find(kv.first);
        if (it != st.removed_by_pid.end()) removed = it->second;
        bypid << _probes[ic].group << "," << _probes[ic].name << "," << kv.first << ","
              << kv.second << "," << removed << "," << frac(removed, kv.second) << "\n";
      }
    }
    bypid.close();

    std::cout << "  wrote " << summary_path << "\n";
    std::cout << "  wrote " << pid_path << "\n";
  }

  long n_events()                const { return _n_events; }
  long n_events_empty_baseline() const { return _n_events_empty_baseline; }

 private:
  struct ProbeStats {
    long n_events_empty = 0;  // events this cut alone left with zero particles
    long n_seen         = 0;  // particles the cut actually judged
    long n_removed      = 0;  // of those, how many it rejected
    std::map<int,long> seen_by_pid;
    std::map<int,long> removed_by_pid;
  };

  static double frac(long a, long b) { return b > 0 ? (double)a / (double)b : 0.0; }

  std::vector<CutProbe>   _probes;
  std::vector<ProbeStats> _stats;
  long _n_events                = 0;
  long _n_particles             = 0;
  long _n_events_empty_nocuts   = 0;
  long _n_events_empty_baseline = 0;
};

#endif
