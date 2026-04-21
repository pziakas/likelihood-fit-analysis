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

    LikelihoodFit fit;

    json input = fit.LoadJSON();

    fit.load_data(fit,input);

    fit.visualize_data(fit,"distribution",c,input);

    fit.perform_fit(fit,c,input);
    
    return EXIT_SUCCESS;
}