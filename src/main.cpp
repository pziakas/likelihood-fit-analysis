// --- C++ includes ---

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <functional>
#include <tuple>
#include <fstream>
#include "json.hpp"
#include "likelihood_fit.h"

// --- ROOT includes ---

#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
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
#include "TSystem.h"

using json = nlohmann::json;

int main(int argc, char **argv)
{
	// --- TApplication is used to show  the plots as in a ROOT script ---

	TApplication app("app", &argc, argv);

	// --- Defining the canvas ---

	TCanvas *c = new TCanvas("c","Plots",800,800);

	TString var;

	// --- Choosing the variable to plot ---

	var_choice(var);

	// --- Choosing the control region --- 

	region_choice(region);

	// --- Opening the json file ---

	std::ifstream in_file("input.json");
	
	json input;

	in_file >> input;

	// --- Opening the file to get the histograms ---
	// --- Remember to put your own path in the JSON file ---

	TString filename(input[region]["filename"].get<std::string>());
	TString path(input["path"].get<std::string>());
	
	TFile *file = TFile::Open(path + filename + var + ".root","read");

	// --- Getting the histograms from the file ---

	TString histo_name(input[region]["histo_name"].get<std::string>());

	get_hists("data_" + histo_name,h_data,file,0);
	get_hists("ZZ_total_" + histo_name,h_ZZ,file,kGreen-9);
	get_hists("WZ_qcd_364253_" + histo_name,h_WZ_qcd,file,40);
	get_hists("WZ_ew_total_" + histo_name,h_WZ_ew,file,kYellow-9);
	get_hists("ttbarV+VVV_total_" + histo_name,h_ttbarV_VVV,file,kMagenta-7);
	get_hists("Fakes_" + histo_name,h_fakes,file,kRed-3);

	histos = {h_ZZ,h_WZ_qcd,h_WZ_ew,h_ttbarV_VVV,h_fakes};

	// --- Getting the singal histo and removing it from the vector ---
	// --- After that the remaining histograms are backgrounds ---

	find_signal_hist(histos);

	// --- Defining useful maps and vectors for the stack and the legend ---


	std::vector <TH1F*> stack_hists_zz = {h_WZ_ew, h_fakes, h_WZ_qcd, h_ttbarV_VVV, h_ZZ};
	std::vector <TH1F*> leg_hists_zz = {h_ZZ, h_ttbarV_VVV, h_WZ_qcd, h_fakes, h_WZ_ew, h_data};

	// --- The wzqcd case ---

	std::vector <TH1F*> stack_hists_wzqcd = {h_WZ_ew, h_ttbarV_VVV, h_ZZ, h_fakes, h_WZ_qcd};
	std::vector <TH1F*> leg_hists_wzqcd = {h_WZ_qcd, h_fakes, h_ZZ, h_ttbarV_VVV, h_WZ_ew, h_data};

	// --- Maps the region to the correct vectors ---

	std::map<std::string,std::vector<std::vector<TH1F*> > > hist_vecs = { { "ZZ", {stack_hists_zz,leg_hists_zz} }, { "WZ_qcd", {stack_hists_wzqcd,leg_hists_wzqcd} } };

	// --- Maps labels to histograms ---

	std::map <TH1F*, std::vector <TString>> label_opts = { {h_ZZ, {"#bf{ZZ}","f"} }, {h_ttbarV_VVV, {"#bf{t#bar{t}V+VVV}","f"} }, {h_WZ_qcd, {"#bf{WZ_{QCD}}","f"} }, {h_fakes, {"#bf{Fake/Non prompt}","f"} }, 

															{h_WZ_ew, {"#bf{WZ_{EW}}","f"} }, {h_data, {"#bf{Data}","ep"} } };



	// --- Drawing the ratio plot ---

	draw_ratio_plot(var,c,hist_vecs[region].at(0),stack,h_data,hist_vecs[region].at(1),label_opts,region,input);


	// --- Drawing the canvas ---

	c->Update();
	c->Draw();


	// --- Asking the user input for printing the canvas ---

	
	print_canvas(c,var + "_distribution");

	
	// --- Changing the size of the canvas and the window ---

	c->SetCanvasSize(700, 500);
	c->GetCanvasImp()->SetWindowSize(700, 500);

	// --- Number of parameters ---

    const int npar =1;
    
    TMinuit *minimizer = new TMinuit(npar);

    // --- Log Likelihood estimation ---

    fit_flag = 0;

    std::vector <double> par = {input[region]["par_guess"].get<double>()};
    std::vector <double> stepSize = {input[region]["par_step"].get<double>()};
    std::vector <double> minVal = {input[region]["par_min"].get<double>()};
    std::vector <double> maxVal = {input[region]["par_max"].get<double>()};
    std::vector <std::string> parName = {input[region]["par_name"].get<std::string>()};


   	best_fit = get_best_fit_values(minimizer,npar,par,stepSize,minVal,maxVal,parName);


    // --- Log Likelihood ratio estimation ---

    fit_flag = 1;

    TMinuit *minimizer2 = new TMinuit(npar);

    best_fit = get_best_fit_values(minimizer2,npar,par,stepSize,minVal,maxVal,parName);

   
    // --- Getting the Log Likelihood ratio ---

    TGraph *ratio = get_likelihood_ratio_plot(minimizer2);


	// --- Estimation of best fit values for the parameter of interest and the nuisance parameter ---

	const int npars = 2;

	// --- Profile likelihood ratio estimation ---

	fit_flag = 2;

	TMinuit *minimizer3 = new TMinuit(npars);

	// --- Adding the information for the second parameter ---

	par.push_back(1.0005);
	stepSize.push_back(0.1);
	minVal.push_back(0.9);
	maxVal.push_back(1.1);
	parName.push_back("x");

	best_fit = get_best_fit_values(minimizer3,npars,par,stepSize,minVal,maxVal,parName);

	
	// --- Getting the minimum and maximum values for the profile likehood ratio plot ---

	double minm = input[region]["mu_min"].get<double>();
	double maxm = input[region]["mu_max"].get<double>();

	
	// --- Getting the profile likelihood ratio ---

	TGraph *profile_ratio = get_profile_likelihood_ratio_plot(minimizer3,minm,maxm,npars,par,stepSize,minVal,maxVal,parName);

	// --- Drawing everything together ---

	draw_ratios(ratio,profile_ratio,c,var,region,input);

	// --- Draws the canvas ---

    c->Update();
    c->Draw();

    // --- Asking for user input to print the canvas ---

	print_canvas(c,var + "_likelihood_plot");

	// --- Done ---

	return 0;

}