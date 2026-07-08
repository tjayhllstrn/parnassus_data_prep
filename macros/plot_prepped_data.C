
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TH1D.h"
#include "TColor.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TROOT.h"
#include "TPie.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"

#include "../src/style.C"

namespace {

std::pair<double, double> get_branch_range(TTree* tree, const std::string& branch_name) {
    double minimum = tree->GetMinimum(branch_name.c_str());
    double maximum = tree->GetMaximum(branch_name.c_str());

    if (!std::isfinite(minimum) || !std::isfinite(maximum)) {
        minimum = 0.0;
        maximum = 1.0;
    }

    if (minimum == maximum) {
        const double pad = (minimum == 0.0) ? 1.0 : std::abs(minimum) * 0.1;
        minimum -= pad;
        maximum += pad;
    }

    const double span = maximum - minimum;
    return {minimum - 0.05 * span, maximum + 0.05 * span};
}

std::unique_ptr<TH1D> make_histogram(TTree* tree,
                                     const std::string& branch_name,
                                     const std::string& hist_name,
                                     int nbins,
                                     double xmin,
                                     double xmax) {
    auto hist = std::make_unique<TH1D>(hist_name.c_str(), "", nbins, xmin, xmax);
    // hist->SetDirectory(nullptr);

    const TString draw_expr = TString::Format("%s>>%s", branch_name.c_str(), hist_name.c_str());
    tree->Draw(draw_expr, "", "goff");
    return hist;
}

void style_overlay_hist(TH1D* hist, int color, int marker_style) {
    hist->SetLineColor(color);
    hist->SetFillColorAlpha(color, 0.7);
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(0.7);
    hist->SetLineWidth(2);
}

void normalize_hist(TH1D* hist) {
    const double integral = hist->Integral("width");
    if (integral > 0.0) {
        hist->Scale(1.0 / integral);
    }
}

void save_plot_objects(const std::string& output_base,
                       TH1D* reco_hist,
                       TH1D* truth_hist,
                       TLegend* legend = nullptr) {
    TFile output_file((output_base + ".root").c_str(), "RECREATE");
    if (!output_file.IsZombie()) {
        if (reco_hist) reco_hist->Write("reco_hist");
        if (truth_hist) truth_hist->Write("truth_hist");
        if (legend) legend->Write("legend");
    }
}

void hsl_to_rgb(double hue, double saturation, double lightness, double& red, double& green, double& blue) {
    const auto hue_to_rgb = [](double p, double q, double t) {
        if (t < 0.0) {
            t += 1.0;
        }
        if (t > 1.0) {
            t -= 1.0;
        }
        if (t < 1.0 / 6.0) {
            return p + (q - p) * 6.0 * t;
        }
        if (t < 1.0 / 2.0) {
            return q;
        }
        if (t < 2.0 / 3.0) {
            return p + (q - p) * (2.0 / 3.0 - t) * 6.0;
        }
        return p;
    };

    if (saturation <= 0.0) {
        red = lightness;
        green = lightness;
        blue = lightness;
        return;
    }

    const double q = lightness < 0.5 ? lightness * (1.0 + saturation)
                                     : lightness + saturation - lightness * saturation;
    const double p = 2.0 * lightness - q;
    red = hue_to_rgb(p, q, hue + 1.0 / 3.0);
    green = hue_to_rgb(p, q, hue);
    blue = hue_to_rgb(p, q, hue - 1.0 / 3.0);
}

int get_pie_color_for_label(const std::string& label) {
    static std::map<std::string, int> color_by_label;
    static std::vector<std::unique_ptr<TColor>> generated_colors;

    const auto found = color_by_label.find(label);
    if (found != color_by_label.end()) {
        return found->second;
    }

    const double golden_ratio = 0.6180339887498949;
    const double hue = std::fmod((color_by_label.size() * golden_ratio), 1.0);
    const double saturation = 0.65;
    const double lightness = 0.55;

    double red = 0.0;
    double green = 0.0;
    double blue = 0.0;
    hsl_to_rgb(hue, saturation, lightness, red, green, blue);

    const int color_index = TColor::GetFreeColorIndex();
    generated_colors.emplace_back(std::make_unique<TColor>(color_index, red, green, blue));
    color_by_label.emplace(label, color_index);
    return color_index;
}

void draw_comparison_plot(TTree* tree,
                          const std::string& reco_branch,
                          const std::string& truth_branch,
                          const std::string& title,
                          const std::string& x_label,
                          const std::string& output_base,
                          double xmin = -999,
                          double xmax = -999,
                          int nbins = 100) {
    if (xmin == -999 || xmax == -999) {
        const auto reco_range = get_branch_range(tree, reco_branch);
        const auto truth_range = get_branch_range(tree, truth_branch);
        const double xmin = std::min(reco_range.first, truth_range.first);
        const double xmax = std::max(reco_range.second, truth_range.second);
    }


    auto reco_hist = make_histogram(tree, reco_branch, "hreco_" + reco_branch, nbins, xmin, xmax);
    auto truth_hist = make_histogram(tree, truth_branch, "htruth_" + truth_branch, nbins, xmin, xmax);

    style_overlay_hist(reco_hist.get(), BLUE, 20);
    style_overlay_hist(truth_hist.get(), MAGENTA, 24);
    normalize_hist(reco_hist.get());
    normalize_hist(truth_hist.get());

    TCanvas canvas(("c_" + output_base).c_str(), title.c_str(), 1200, 900);
    canvas.cd();

    reco_hist->SetTitle(title.c_str());
    reco_hist->GetXaxis()->SetTitle(x_label.c_str());
    reco_hist->GetYaxis()->SetTitle("Normalized entries");
    reco_hist->Draw("hist");
    truth_hist->Draw("hist same");

    TLegend legend(0.68, 0.78, 0.90, 0.90);
    legend.AddEntry(reco_hist.get(), (reco_branch).c_str(),"f");
    legend.AddEntry(truth_hist.get(), (truth_branch).c_str(),"f");
    legend.Draw();

    canvas.SaveAs((output_base + ".png").c_str());
    save_plot_objects(output_base, reco_hist.get(), truth_hist.get(), &legend);
    // canvas.SaveAs((output_base + ".pdf").c_str());
}

void label_integer_bins(TH1D* hist, const std::map<int, long long>& counts) {
    int bin = 1;
    for (const auto& entry : counts) {
        hist->GetXaxis()->SetBinLabel(bin, TString::Format("%d", entry.first));
        hist->SetBinContent(bin, entry.second);
        ++bin;
    }
}

std::unique_ptr<TH1D> build_category_histogram(TTree* tree,
                                               const std::string& branch_name,
                                               const std::string& hist_name) {
    std::map<int, long long> counts;
    TTreeReader reader(tree);
    TTreeReaderArray<int> values(reader, branch_name.c_str());

    while (reader.Next()) {
        for (const auto value : values) {
            counts[value]++;
        }
    }

    const int nbins = std::max(1, static_cast<int>(counts.size()));
    auto hist = std::make_unique<TH1D>(hist_name.c_str(), "", nbins, 0, nbins);
    hist->SetDirectory(nullptr);
    hist->GetXaxis()->SetTitle(branch_name.c_str());
    hist->GetYaxis()->SetTitle("Entries");
    label_integer_bins(hist.get(), counts);
    return hist;
}

auto draw_pie_chart(TH1D* hist,
                    const std::string& title,
                    const std::string& output_base) {
    std::vector<double> values;
    std::vector<int> colors;
    std::vector<std::string> labels_storage;
    std::vector<const char*> labels;

    const int nbins = hist->GetNbinsX();
    values.reserve(nbins);
    colors.reserve(nbins);
    labels_storage.reserve(nbins);
    labels.reserve(nbins);

    for (int bin = 1; bin <= nbins; ++bin) {
        const std::string label = hist->GetXaxis()->GetBinLabel(bin);
        values.push_back(hist->GetBinContent(bin));
        colors.push_back(get_pie_color_for_label(label));
        labels_storage.push_back(label);
        labels.push_back(labels_storage.back().c_str());
    }


    auto pie = std::make_unique<TPie>(
        ("pie_" + output_base).c_str(),
        title.c_str(),
        nbins,
        values.data(),
        colors.data(),
        labels.data());
    pie->SetLabelFormat("%txt (%perc)");
    pie->SetTextFont(42);
    pie->SetTextSize(0.03);
    pie->SetTextColor(kBlack);
    pie->SetCircle(0.5, 0.5, 0.35);
    for (int slice = 0; slice < nbins; ++slice) {
        pie->SetEntryLineColor(slice, kWhite);
        pie->SetEntryLineWidth(slice, 1);
    }
    // pie->Draw("SC");

    return pie;
}

void draw_pie_comparison(TH1D* reco_hist,
                         TH1D* truth_hist,
                         const std::string& title,
                         const std::string& output_base) {
    TCanvas canvas(("c_" + output_base).c_str(), title.c_str(), 1400, 700);
    canvas.Divide(2, 1);

    canvas.cd(1);
    gPad->SetMargin(0.05, 0.05, 0.10, 0.10);
    auto pie_reco = draw_pie_chart(reco_hist, title + " reco", output_base + "_reco");
    pie_reco->Draw("SC");

    canvas.cd(2);
    gPad->SetMargin(0.05, 0.05, 0.10, 0.10);
    auto pie_truth = draw_pie_chart(truth_hist, title + " truth", output_base + "_truth");
    pie_truth->Draw("SC");

    canvas.SaveAs((output_base + ".png").c_str());
    save_plot_objects(output_base, reco_hist, truth_hist);
    // canvas.SaveAs((output_base + ".pdf").c_str());
}

} // namespace

void plot_prepped_data(const char* filename){
    SetReleaseNoteStyle();
    DefineColors();

    std::vector<std::string> vars_to_plot = {"pt", "phi", "eta", "vx", "vy", "vz"};
    std::vector<std::string> x_labels = {"p_{T} [GeV/c]", "#phi [rad]", "#eta", "v_{x} [cm]", "v_{y} [cm]", "v_{z} [cm]"};
    std::vector<std::pair<double, double>> val_ranges {{-0.1,2}, {-3.14159, 3.14159}, {-1, 5}, {-4, 4}, {-4, 4}, {-10, 10}};

    std::unique_ptr<TFile> input_file(TFile::Open(filename, "READ"));
    if (!input_file || input_file->IsZombie()) {
        Error("plot_prepped_data", "Could not open input file: %s", filename);
        return;
    }

    TTree* tree = nullptr;
    input_file->GetObject("evt_tree", tree);
    if (!tree) {
        Error("plot_prepped_data", "Could not find tree evt_tree in file: %s", filename);
        return;
    }

    std::string outdir = "plots/" + std::string(filename);
    gSystem->mkdir(outdir.c_str(), true);

    for (size_t i = 0; i < vars_to_plot.size(); ++i) {    
        draw_comparison_plot(tree,
                            "reco_" + vars_to_plot[i],
                            "truth_" + vars_to_plot[i],
                            "reco vs truth " + vars_to_plot[i],
                            x_labels[i],
                            outdir + "/compare_" + vars_to_plot[i],
                            val_ranges[i].first,
                            val_ranges[i].second);
    }

    draw_comparison_plot(tree,
                         "nreco",
                         "ntruth",
                         "nreco vs ntruth",
                         "multiplicity",
                         outdir + "/compare_nreco_ntruth",
                         -999,-999,20);

    auto reco_class_hist = build_category_histogram(tree, "reco_class", "hreco_class");
    auto truth_class_hist = build_category_histogram(tree, "truth_class", "htruth_class");
    auto reco_pid_hist = build_category_histogram(tree, "reco_pid", "hreco_pid");
    auto truth_pid_hist = build_category_histogram(tree, "truth_pid", "htruth_pid");

    draw_pie_comparison(reco_class_hist.get(),
                        truth_class_hist.get(),
                        "classID",
                        outdir + "/classID_pie_comparison");
    draw_pie_comparison(reco_pid_hist.get(),
                        truth_pid_hist.get(),
                        "PID",
                        outdir + "/PID_pie_comparison");
}