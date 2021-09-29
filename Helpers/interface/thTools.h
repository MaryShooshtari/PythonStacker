#ifndef TH_TOOLS_H
#define TH_TOOLS_H

#include <TROOT.h>
#include <TH1.h>
#include <THStack.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TColor.h>

void normalizeHistograms(std::vector<TH1D*>& histos);

#endif