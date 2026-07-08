#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"

TStyle* ReleaseNoteStyle() {

    TStyle *s = new TStyle("ReleaseNoteStyle","Poster plotting style");

    // Canvas
    s->SetCanvasColor(0);
    s->SetCanvasBorderMode(0);
    s->SetCanvasDefW(1200);
    s->SetCanvasDefH(1000);

    // Pads
    s->SetPadColor(0);
    s->SetPadBorderMode(0);
    s->SetPadTopMargin(0.07);
    s->SetPadBottomMargin(0.12);
    s->SetPadLeftMargin(0.13);
    s->SetPadRightMargin(0.05);

    // Frame
    s->SetFrameBorderMode(0);
    s->SetFrameLineWidth(0);

    // Fonts
    s->SetTextFont(42);
    s->SetLabelFont(42,"XYZ");
    s->SetTitleFont(42,"XYZ");

    //Title
    s->SetTitleFillColor(0);    // removes grey background
    s->SetTitleBorderSize(0);

    // Axis sizes
    s->SetTitleSize(0.05,"XYZ");
    s->SetLabelSize(0.045,"XYZ");

    // Axis offsets
    s->SetTitleOffset(1.1,"X");
    s->SetTitleOffset(1.3,"Y");

    //Set max digits on y axis
    //    s->SetMaxDigits(4,"Y");

    // Ticks
    s->SetPadTickX(0);
    s->SetPadTickY(0);

    // Error bars
    s->SetErrorX(0);

    // Histogram defaults
    s->SetHistLineWidth(2);

    // Turn off stats box
    s->SetOptStat(0);

    // Legend style
    s->SetLegendBorderSize(0);

    // Main Title
    s->SetTitleBorderSize(0);
    s->SetTitleX(0.5);
    s->SetTitleAlign(23);

    return s;
}

void SetReleaseNoteStyle() {
    TStyle* s = ReleaseNoteStyle();
    gROOT->SetStyle("ReleaseNoteStyle");
    gROOT->ForceStyle();

    TGaxis::SetMaxDigits(3);
}

int BLUE;
int MAGENTA;
void DefineColors() {
    BLUE = TColor::GetColor(68, 119, 170);
    MAGENTA = TColor::GetColor(170, 51, 119);
}
