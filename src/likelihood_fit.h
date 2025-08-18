#ifndef LIKELIHOOD_FIT_H
#define LIKELIHOOD_FIT_H


// --- C++ includes ---

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <functional>
#include <fstream>
#include "json.hpp"
#include <tuple>

// --- ROOT includes ---

#include "TH1.h"
#include "TFile.h"
#include "TF1.h"
#include "TString.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TApplication.h"
#include "TLine.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Minuit2/MnUserParameterState.h"
#include "Math/IFunction.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include "TGraph.h"
#include "TCanvasImp.h"

using json = nlohmann::json;

// --- Histograms ---

extern TH1F *h_data;
extern TH1F *h_ZZ;
extern TH1F *h_WZ_qcd;
extern TH1F *h_WZ_ew;
extern TH1F *h_fakes;
extern TH1F *h_ttbarV_VVV;
extern TH1F *h_signal;


// --- Region ---

extern std::string region;

// --- Flags ---

extern int fit_flag;

// --- Global vectors  ---

extern std::vector <double> best_fit;
extern std::vector <TH1F*> histos;


// --- Functions ---

void draw_line(int color, int line_width, double xmin, double ymin, double xmax, double ymax);
void draw_logo(double xmin, double xmax1, double xmax2, double xmax3,std::string region, const json &input);
void draw_legend(const std::vector <TH1F*> &hists, const std::map <TH1F*,std::vector<TString>> &leg_map);
void var_choice(TString &var);
void get_hists(const TString &name, TH1F *&h, TFile *file, int fill_color);
void style_hists(TH1F *&h, int fill_color);
void draw_ratio_plot(const TString &var, TCanvas *c, const std::vector <TH1F*> &hists, THStack *&s, TH1F *h_data,const std::vector <TH1F*> &hists_leg, const std::map <TH1F*,std::vector<TString>> &leg_map,std::string region, const json &input);
double log_likelihood(TH1F *h_bkg1, TH1F *h_bkg2, TH1F *h_bkg3, TH1F *h_bkg4, TH1F *h_signal, TH1F *h_data, int bin, const double *par, int flag);
void fcn(int &npar, double *deriv, double &f, double *par, TH1F *h_bkg1, TH1F *h_bkg2, TH1F *h_bkg3, TH1F *h_bkg4, TH1F *h_signal, TH1F *h_data, int bin, int flag);
std::vector<double> get_best_fit_values(TMinuit *&minimizer,const int npar, const std::vector <double> &par, const std::vector <double> &stepSize, const std::vector <double> &minVal,const std::vector <double> &maxVal, const std::vector <std::string> &parName);
TGraph* get_likelihood_ratio_plot(TMinuit *&minimizer);
void draw_sigma(double x, double y);
void style_ratios(TGraph *plot, int color, std::string region, const json &input);
TGraph* get_profile_likelihood_ratio_plot(TMinuit *&minimizer, double minm, double maxm, const int npars, std::vector <double> &par, const std::vector <double> &stepSize, const std::vector <double> &minVal,const std::vector <double> &maxVal, const std::vector <std::string> &parName);
std::tuple<double,double,double> quadratic_fit(TGraph *plot, std::string region, const json &input);
void draw_ratios(TGraph *ratio, TGraph *profile_ratio, TCanvas *canvas, const TString &var, std::string region, const json &input);
void print_canvas(TCanvas *c, TString image);
void draw_likelihood_legend(TGraph *ratio, TGraph *profile_ratio);
void region_choice(std::string &region);
void find_signal_hist(std::vector<TH1F*> &histos);



#endif