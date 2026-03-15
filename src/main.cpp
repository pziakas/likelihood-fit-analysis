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

	auto c = std::make_unique<TCanvas>("c", "Canvas", 800, 800);

	std::ifstream in_file("input.json");
	
	json input;

	in_file >> input;

    LikelihoodFit fit;

    fit.ShowVars(input);

    try
    {
        fit.SetVar(input);

        fit.SetRegion(input);

        fit.GetVectors(input);

        fit.FileOpen(input);

        fit.GetHistos(input);
    }

    catch(const std::exception& e)
    {
        std::cerr << "ERROR: " <<  e.what() << std::endl;
        return EXIT_FAILURE;
    }

    fit.StyleHistos();


    fit.draw_ratio_plot(c,input);


    fit.print_image(c,"distribution",input);

    // --- Changing the size of the canvas and the window for the likelihood plot ---

	c->SetCanvasSize(700, 500);
	c->GetCanvasImp()->SetWindowSize(700, 500);

    // --- Log Likelihood estimation ---

    fit.SetFitFlag(0);


    fit.LoadParameters(input);

    fit.SetNpar(1);

    LikelihoodFit::SetInstance(fit);

    fit.SetBestFitValues();

    // --- Log Likelihood Ratio estimation ---

    fit.SetFitFlag(1);

    // --- Getting the Log Likelihood Ratio ---

    auto ratio = std::unique_ptr<TGraph>();

    try
    {
       ratio = fit.get_likelihood_ratio_plot(fit);
    } 
    
    catch(const std::exception& e)
    {
        std::cerr << "ERROR: " <<  e.what() << std::endl;
        return EXIT_FAILURE;
    }
   
    fit.SetNpar(2);

    // --- Profile Likelihood Ratio estimation ---

    fit.SetFitFlag(2);

    // --- Adding the initial values for the second parameter ---

    fit.AddToParameters(1.0005,0.1,0.9,1.1,"x");

    fit.SetBestFitValues();

    // --- Getting the Profile Likelihood Ratio ---

    auto profile_ratio = fit.get_profile_likelihood_ratio_plot(input);

	fit.draw_ratios(ratio,profile_ratio,c,input);

    fit.print_image(c,"likelihood",input);
    
    return EXIT_SUCCESS;
}