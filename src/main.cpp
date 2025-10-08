// --- C++ includes ---

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <functional>
#include <tuple>
#include <fstream>
#include <memory>
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

    // --- Defining a canvas ---

	auto c = std::make_unique<TCanvas>("c", "Canvas", 800, 800);

    // --- Creating a LikelihoodFit object ---

    LikelihoodFit fit;

    // --- Setting the variable to be plotted ---

    fit.SetVar();

    // --- Choosing the region that will be plotted ---

    fit.SetRegion();

    // --- Opening the JSON file ---

	std::ifstream in_file("input.json");
	
	json input;

	in_file >> input;

    // --- Loading the usefull vectors ---

    fit.GetVectors(input);

    // --- Getting the path from the JSON file ---

    fit.FileOpen(input);

    // --- Getting the histograms from the file ---

    fit.GetHistos(input);

    // --- Checking if the histograms are correctly loaded ---

    fit.CheckHistos();

    // --- Styling the histograms ---

    fit.StyleHistos();

    // --- Drawing the ratio plot ---

    fit.draw_ratio_plot(c,input);

    // --- Drawing the canvas ---

    c->Update();
    c->Draw();

    // --- Asking the user if they want to prin the stack ---

    fit.print_image(c,"distribution");

    // --- Changing the size of the canvas and the window ---

	c->SetCanvasSize(700, 500);
	c->GetCanvasImp()->SetWindowSize(700, 500);

    // --- Log Likelihood estimation ---

    fit.SetFitFlag(0);

    // --- Loading the parameters ---

    fit.LoadParameters(input);

    // --- Setting the number of parameters for the fit ---

    fit.SetNpar(1);

    // --- Setting the instant ---

    LikelihoodFit::SetInstance(fit);

    // --- Performing the minimization ---

    fit.SetBestFitValues();

    // --- Log Likelihood Ratio estimation ---

    fit.SetFitFlag(1);

    // --- Getting the Log Likelihood Ratio ---

    TGraph *ratio = fit.get_likelihood_ratio_plot(fit);

    // --- Setting the number of parameters for the fit ---

    fit.SetNpar(2);

    // --- Profile Likelihood Ratio estimation ---

    fit.SetFitFlag(2);

    // --- Adding the initial values for the second parameter ---

    fit.AddToParameters(1.0005,0.1,0.9,1.1,"x");

    // --- Performing the minimization to obtain best fit values for both parameters ---
    
    fit.SetBestFitValues();

    // --- Getting the Profile Likelihood Ratio ---

    TGraph *profile_ratio = fit.get_profile_likelihood_ratio_plot(input);

    // --- Drawing everything together ---

	fit.draw_ratios(ratio,profile_ratio,c,input);

	// --- Draws the canvas ---

    c->Update();
    c->Draw();

    // --- Asking the user if they want to prin the likelihood plot ---

    fit.print_image(c,"likelihood");

    // --- Done ---
    
    return 0;
}