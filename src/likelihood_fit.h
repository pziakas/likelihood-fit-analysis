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
#include <memory>

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

class LikelihoodFit
{
    // --- Defining some usefull vectors as members of the class ---

    std::vector <TH1F*> histos;
    std::vector <TString> histo_names;
    std::vector <TString> leg_names;
    std::vector <int> fill_color;

    // --- Some usefull TStrings ---

    TString var;
    TString region;

    // --- File to be opened ---

    TFile *file;

    // --- Defining some usefull integers ---

    int fit_flag;
    int npar;

    // --- Defining some helpful vectors for the parameters of the fit ---

    std::vector<double> par;
    std::vector<double> stepSize;
    std::vector<double> minVal;
    std::vector<double> maxVal;
    std::vector<std::string> parName;
    std::vector<double> best_fit;

    // --- Static instance member ---

    static LikelihoodFit* instance;
    
    public:

    LikelihoodFit();
    ~LikelihoodFit();
    void FileOpen(const json &input);
    void SetRegion();
    TString GetRegion();
    void SetVar();
    void GetHistos(const json &input);
    void StyleHistos();
    void CheckHistos();
    void draw_logo(double xmin, double xmax1, double xmax2, double xmax3, const json &input);
    void draw_legend();
    void draw_line(int color, int line_width, double xmin, double ymin, double xmax, double ymax);
    void draw_ratio_plot(std::unique_ptr<TCanvas> &c, const json &input);
    void GetVectors(const json &input);
    std::vector <int> vecs_from_json_int(const json &input, const std::string &name);
    std::vector <TString> vecs_from_json_string(const json &input, const std::string &name);
    void print_image(std::unique_ptr<TCanvas> &c, const TString &mode);
    void SetFitFlag(const int &fit_flag);
    double log_likelihood(int bin, double param[], int flag);
    void LoadParameters(const json &input);
    static void fcn(int &npar, double *deriv, double &f, double *param, int flag);
    std::vector<double> get_best_fit_values(TMinuit *minimizer);
    void SetNpar(int npar);
    int GetNpar();
    static void SetInstance(LikelihoodFit &fit);
    TGraph* get_likelihood_ratio_plot(LikelihoodFit &fit);
    void AddToParameters(double par1, double par2, double par3, double par4, std::string par5);
    TGraph* get_profile_likelihood_ratio_plot(const json &input);
    TMinuit* SetBestFitValues();
    void style_ratios(TGraph *plot, int color, const json &input);
    void draw_sigma(double x, double y);
    std::tuple<double,double,double> quadratic_fit(TGraph *plot, const json &input);
    void draw_likelihood_legend(TGraph *ratio, TGraph *profile_ratio);
    void draw_ratios(TGraph *ratio, TGraph *profile_ratio,std::unique_ptr<TCanvas> &c, const json &input);

};



#endif