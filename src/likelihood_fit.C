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
#include "likelihood_fit.h"

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

// ########################
// #                      #
// # Function definitions #
// #                      #
// ########################

// --- Static pointer of the class ---

LikelihoodFit* LikelihoodFit::instance = nullptr;

// --- Setter to set the instance ---

//************************************************************************
void LikelihoodFit::SetInstance(LikelihoodFit &fit)
//************************************************************************
{
    instance = &fit;
}

// --- LikelihoodFit Constructor ---

//************************************************************************
LikelihoodFit::LikelihoodFit()
//************************************************************************
{
    std::cout << "You have created a LikelihoodFit object!" << std::endl << std::endl;
}

// --- LikelihoodFit Destructor ---

//************************************************************************
LikelihoodFit::~LikelihoodFit()
//************************************************************************
{
    for(auto &histo : histos)
    {
        delete histo;
    }
    
    histos.clear();

    file->Close();
    delete file;

    std::cout << "You have deleted a LikelihoodFit object" << std::endl;
}

// --- This function is opening the file and checks if it has opened successfully ---

//************************************************************************
void LikelihoodFit::FileOpen(const json &input)
//************************************************************************
{
    TString filename(input[region]["filename"].get<std::string>());
	TString path(input["path"].get<std::string>());
    
    file = TFile::Open(path + filename + var + ".root","read");

    if(!file)
    {
        std::cout << std::endl << "Could not open the file! Exiting..." << std::endl;
        std::exit(1);
    }

    std::cout << std::endl << "File " << file->GetName() << " opened successfully!" << std::endl;
}

// --- This is a function to set the region ---

//************************************************************************
void LikelihoodFit::SetRegion()
//************************************************************************
{
    TString region;

    std::cout << "Enter the region you want to plot: ";
    std::cin >> region;

    while(region != "ZZ" && region != "WZ_qcd")
    {
        std::cout << "Not valid region!" << std::endl;
        std::cout << "You can choose between: ZZ and WZ_qcd!" << std::endl;
        std::cout << "Enter the region you want to plot: ";
        std::cin >> region;
    }

    this->region = region;
}

// --- This is a getter function for the region ---

//************************************************************************
TString LikelihoodFit::GetRegion()
//************************************************************************
{
    return region;
}

// --- This is a function to set the variable to be plotted ---

//************************************************************************
void LikelihoodFit::SetVar()
//************************************************************************
{
    TString var;

    std::cout << "Enter the variable you want to plot: ";
    std::cin >> var;

    while(var != "mwz" && var != "sum3pt" && var != "mtwz")
	{
		std::cout << "Not valid variable!" << std::endl;
        std::cout << "You can choose between: mwz, sum3pt or mtwz!" << std::endl;
		std::cout << "Enter the variable you want to plot: ";
		std::cin >> var;
	}

   this->var = var;
}

// --- This is a function that gets all the histograms --

//************************************************************************
void LikelihoodFit::GetHistos(const json &input)
//************************************************************************
{
    TString histo_label(input[region]["histo_label"].get<std::string>());
    
    for(auto &histo_name : histo_names)
    {
        TString name = histo_name + histo_label;
        
        // --- Getting the histogram from the file ---
        
        TH1F *histo = dynamic_cast<TH1F*>(file->Get(name));

        // --- Setting the name of the histogram ---
        
        histo->SetName(histo_name);

        // --- Filling the histogram vector ---
        
        histos.push_back(histo);
    }
}

// --- This function styles the histograms ---

//************************************************************************
void LikelihoodFit::StyleHistos()
//************************************************************************
{
    // --- Looping over all the histograms ---
    
    for(int i = 0; i < histos.size(); i++)
    {
        TString name = histos.at(i)->GetName();

        // --- Styling the data histogram ---

        if(name.Contains("data"))
        {
            histos.at(i)->SetMarkerStyle(20);
            histos.at(i)->SetMarkerColor(kBlack);
            histos.at(i)->SetLineColor(kBlack);
        }

        // --- All the predefined vectors have the same size ---

        histos.at(i)->SetFillColor(fill_color.at(i));
     
    }

}

// --- This function checks if the histograms have been correctly taken from the file ---

//************************************************************************
void LikelihoodFit::CheckHistos()
//************************************************************************
{
    for(const auto &histo : histos)
    {
        if(!histo)
        {
            std::cout << std::endl << "Histogram could not be successfully opened! Exiting..." << std::endl;
            std::exit(1);
        }

        std::cout << std::endl << "Histogram " << histo->GetName() << " opened successfully!" << std::endl;
    }
}

// --- This function draws the ATLAS logo ---

//************************************************************************
void LikelihoodFit::draw_logo(double x, double y1, double y2, double y3, const json &input)
//************************************************************************
{
    TString region_logo(input[region]["region"].get<std::string>());
    
    TLatex *t = new TLatex(x,y1,"#it{ATLAS} #bf{work in progress}");
    TLatex *t1=new TLatex(x,y2,"#bf{#sqrt{s}=13 TeV, #intLdt = 139 fb^{-1}}");
    TLatex *t2=new TLatex(x,y3,region_logo);
    
    t->SetTextSize(0.04);
    t1->SetTextSize(0.035);
    t2->SetTextSize(0.035);
    
    t->Draw();
    t1->Draw();
    t2->Draw();
}

// --- This function draws the legend for the stack ---

//************************************************************************
void LikelihoodFit::draw_legend()
//************************************************************************
{
	TLegend *leg = new TLegend(0.6,0.55,0.85,0.85); 
	
    // --- The two vectors have the same size ---   

    for(int i = 0; i < histos.size(); i++)
    {
        if(leg_names.at(i).Contains("Data")) leg->AddEntry(histos.at(i),leg_names.at(i),"ep");
        else leg->AddEntry(histos.at(i),leg_names.at(i),"f");
    }

	leg->SetBorderSize(0);
	leg->Draw();
}

// --- This function plots the line of the ratio plot ---

//************************************************************************
void LikelihoodFit::draw_line(int color, int line_width, double xmin, double ymin, double xmax, double ymax)
//************************************************************************
{
	
	TLine *l = new TLine(xmin,ymin,xmax,ymax);
    l->SetLineColor(color);
	l->SetLineWidth(line_width);
    l->SetLineStyle(9);
    l->Draw("same");
}

// --- This function will draw the ratio plot ---

//************************************************************************
void LikelihoodFit::draw_ratio_plot(std::unique_ptr<TCanvas> &c, const json &input)
//************************************************************************
{
    // --- Defining the stack ---

     THStack *s = new THStack("stack",";; Counts");
	
	// --- Getting the stack limits from the json file ---

	double stack_max = input[region][var]["stack_upper_lim"].get<double>();
    double stack_min = input[region][var]["stack_lower_lim"].get<double>();


	// --- Filling the stack ---
    // --- First histogram is data histosgram, will be plotted later ---

	for(int i = 1; i < histos.size(); i++) s->Add(histos.at(i));

	// --- Drawing the upper pad
	
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.05); 
    pad1->Draw();
    pad1->cd();           

    // --- Drawing the stack ---

    s->SetMinimum(stack_min);
    s->SetMaximum(stack_max);
    s->Draw("hist");

    s->GetXaxis()->SetTitle("");
    s->GetXaxis()->SetLabelOffset(999);
    s->GetXaxis()->SetLabelSize(0);

    // --- No stats on the plots ---

    gStyle->SetOptStat(0); 

    // --- Draw the logo ---

    double logo_y1 = input[region][var]["atlas_logo_y1"].get<double>();
    double logo_y2 = input[region][var]["atlas_logo_y2"].get<double>();
    double logo_y3 = input[region][var]["atlas_logo_y3"].get<double>();
    double logo_x = input[region]["stack_logo_x"].get<double>();

    draw_logo(logo_x,logo_y1,logo_y2,logo_y3,input);

    // --- Draw the legend ---
    
    draw_legend();

    // --- Draw data ---

    histos.at(0)->Draw("e1x0p same");

    // --- Logarithmic scale ---

    pad1->SetLogy();

    // --- Going back to the main canvas and drawing the lower pad ---

    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd(); 

    // --- Clone the histograms for the ratio plot ---

    TH1F *h_Data = dynamic_cast<TH1F*>(histos.at(0)->Clone("h_Data"));
    TH1F *h_sum = dynamic_cast<TH1F*>(histos.at(1)->Clone("h_sum"));

    // --- Not taking again into account the second histo for the background sum ---
    // --- This is because it was already copied above ---
    // --- Otherwise its contribution is counted twice ---

    for(int i = 2; i < histos.size(); i++) h_sum->Add(histos.at(i));

    // --- Divide the histograms and draw the ratio ---

    h_Data->SetTitle("");
    h_Data->SetStats(0);      
    h_Data->Divide(h_sum);
    h_Data->SetMarkerStyle(20);
    h_Data->SetLineColor(kBlack);
    h_Data->Draw("ex0 ep");   

    // --- Ratio plot limits ---

    h_Data->SetMinimum(0.5);
    h_Data->SetMaximum(1.5);

    // --- X axis ---

    h_Data->GetXaxis()->SetTitle(var + " [GeV]");
    h_Data->GetXaxis()->SetTitleSize(20);
    h_Data->GetXaxis()->SetTitleFont(43);
    h_Data->GetXaxis()->SetTitleOffset(1);
    h_Data->GetXaxis()->SetLabelFont(43); 
    h_Data->GetXaxis()->SetLabelSize(15);  

    // --- Y axis ---

    h_Data->GetYaxis()->SetTitle("Data/Prefit");
    h_Data->GetYaxis()->SetNdivisions(505);
    h_Data->GetYaxis()->SetTitleSize(20);
    h_Data->GetYaxis()->SetTitleFont(43);
    h_Data->GetYaxis()->SetTitleOffset(1.55);
    h_Data->GetYaxis()->SetLabelFont(43); 
    h_Data->GetYaxis()->SetLabelSize(15);  

    // --- All the histograms have the same limits in their X axis ---
    // --- The limits can be computed like that ---

    int first_bin = histos.at(0)->FindFirstBinAbove(-10);
    int last_bin = histos.at(0)->FindLastBinAbove(-10);

    double minx = histos.at(0)->GetBinLowEdge(first_bin);
    double maxx = histos.at(0)->GetBinLowEdge(last_bin) + histos.at(0)->GetBinWidth(last_bin);

    draw_line(kBlack,1,minx,1,maxx,1);
 }

// --- This function will turn the JSON file string vectors in to TString vectors ---

//************************************************************************
std::vector <TString> LikelihoodFit::vecs_from_json_string(const json &input, const std::string &name)
//************************************************************************
{
    std::vector <TString> vec;

    for(const auto &elements : input[region][name])
    {
        TString element(elements.get<std::string>());
        vec.push_back(element);
    }

    return vec;
}

// --- This function will turn the JSON file int vectors in to TString vectors ---

//************************************************************************
std::vector <int> LikelihoodFit::vecs_from_json_int(const json &input, const std::string &name)
//************************************************************************
{
    std::vector <int> vec;

    for(const auto &elements : input[region][name])
    {
        int element = elements;
        vec.push_back(element);
    }

    return vec;
}

// --- This function will set all the vectors ---

//************************************************************************
void LikelihoodFit::GetVectors(const json &input)
//************************************************************************
{
    histo_names = vecs_from_json_string(input,"histo_names");
    leg_names = vecs_from_json_string(input,"leg_names");
    fill_color = vecs_from_json_int(input,"fill_colors");
}

// --- This function is asking the user if they want to print the canvas in a pdf form ---

//************************************************************************
void LikelihoodFit::print_image(std::unique_ptr<TCanvas> &c, const TString &mode)
//************************************************************************
{
    std::string ans;
    
    std::cout << std::endl << "Do you want to print the canvas [y/n]: ";
    std::cin >> ans;

    while(ans != "y" && ans!= "n")
    {
        std::cout << "You have entered an invalid choice!" << std::endl;
        std::cout << "Do you want to print the canvas [y/n]: ";
        std::cin >> ans;
    }

    if(ans == "y")
    {
        c->Print(var + "_" + mode + ".png");
    }

}

// --- This function sets the fit flag for the different fit cases ---

//************************************************************************
void LikelihoodFit::SetFitFlag(const int &fit_flag)
//************************************************************************
{
    this->fit_flag = fit_flag;
}

// --- This function sets the number of parameters for the fit ---

//************************************************************************
void LikelihoodFit::SetNpar(int npar)
//************************************************************************
{
    this->npar = npar;
}

// --- This function gives the number of parameters for the fit ---

//************************************************************************
int LikelihoodFit::GetNpar()
//************************************************************************
{
    return npar;
}

// --- This function gives the binned likelihood value ---
// --- The value that is returned is different for each bin ---
// --- The flag corresponds to the calculation with or withour the nuisance parameter ---

//************************************************************************
double LikelihoodFit::log_likelihood(int bin, double param[], int flag)
//************************************************************************
{
    double val;

    // --- Defining the needed fit functions ---
    
    TF1 *poisson = new TF1("p","TMath::Poisson(x,[0])",0,1e5);
    TF1 *gauss = new TF1("g","gaus(0)",-10,10);

    double s = 0;
    double b = 0;

    // --- The histogram that contains the name of region is considered as signal ---
    // --- The rest is considered as background ---
    
    for(int i = 1; i < histos.size(); i++)
    {
        TString name = histos.at(i)->GetName();

        if(name.Contains(region))
        {
            s = histos.at(i)->GetBinContent(bin);
            continue;
        }

        b += histos.at(i)->GetBinContent(bin);
    }

    // --- The luminosity systematic is inserted in the likelihood as a nuisance parameter p[0] which can change the number of total predicted events ---
    // --- This parameter is known to follow a Gaussian distribution, is centered around 1 and has a realative error of 2% ---

    gauss->SetParameter(0,1/(sqrt(TMath::TwoPi())*(0.02)));
    gauss->SetParameter(1,1);
    gauss->SetParameter(2,0.02);

    // --- Corresponds to the non-systematics case log-likelihood estimation ---
 
    if(flag == 0)
    {
        poisson->SetParameter(0,( (param[0]*s) + b ) );
       
        val = log(poisson->Eval(histos.at(0)->GetBinContent(bin)));
    } 

    // --- Corresponds to the calculation with the systematic uncertainty ---

    else if(flag == 1)
    {
        poisson->SetParameter(0,param[1]*( (param[0]*s) + b ) );

        val = log(poisson->Eval(histos.at(0)->GetBinContent(bin))) + log(gauss->Eval(param[1]));
    } 

    return val;

}

// --- This is the function that is going to be given to the minimizer ---
// --- It can't get different arguments, otherwise TMinuit is not working ---

//************************************************************************
void LikelihoodFit::fcn(int &npar, double *deriv, double &f, double *param, int flag)
//************************************************************************
{
    
    double lnL = 0;

    int size = instance->best_fit.size();

    double vals[size];

    for(int i = 0; i < size; i++) vals[i] = instance->best_fit.at(i);

    // --- All of the histograms have the same number of bins ---

    int nbins = instance->histos.at(0)->GetNbinsX();

    for(int i = 1; i <= nbins; i++)
    {
        // --- Binned likelihood fit ---

        if(instance->fit_flag == 0)
        {
            lnL = lnL + instance->log_likelihood(i, param, 0);
        }

        // --- Binned likelihood ratio fit ---

        else if(instance->fit_flag == 1) 
        {
            vals[0] = instance->best_fit.at(0);
            
            lnL = lnL + instance->log_likelihood(i, param, 0) - instance->log_likelihood(i, vals, 0);
        
        }

         // --- Likelihood ratio fit with one nuisance parameter ---

        else if(instance->fit_flag == 2)
        {
            lnL = lnL + instance->log_likelihood(i, param, 1);
        } 
    }

    f = -2*lnL;
}

// --- This function will load the parameters for the fit ---

//************************************************************************
void LikelihoodFit::LoadParameters(const json &input)
//************************************************************************
{
    // --- Loading the vectors ---
    
    par = {input[region]["par_guess"].get<double>()};
    stepSize = {input[region]["par_step"].get<double>()};
    minVal = {input[region]["par_min"].get<double>()};
    maxVal = {input[region]["par_max"].get<double>()};
    parName = {input[region]["par_name"].get<std::string>()};
}

// --- This function gets a TMinuit minimizer of npar parameters and returns a vector of the best fit values ---

//************************************************************************
std::vector<double> LikelihoodFit::get_best_fit_values(TMinuit *minimizer)
//************************************************************************
{
    std::vector<double> values;

    // --- Setting the minimization function ---
    
    minimizer->SetFCN(LikelihoodFit::fcn);

    // --- Setting the parameters ---

    for(int i = 0; i < npar; i++) minimizer->DefineParameter(i,parName[i].c_str(),par[i],stepSize[i],minVal[i],maxVal[i]);

    // --- The algorithm that is called for the minimization ---
    
    minimizer->Migrad();
   
    double outpar[npar];
    double error[npar];

    for(int i=0;i<npar;i++)
    {
       minimizer->GetParameter(i,outpar[i],error[i]);

       std::cout << parName[i] << ": " << outpar[i] << " +/- " << error[i] << std::endl;

       // --- Keeping the best fit values in the vector ---

       values.push_back(outpar[i]);  
    
    }
    
    return values;
   
}

// --- This function plots the likelihood ratio estimation ---

//************************************************************************
TGraph* LikelihoodFit::get_likelihood_ratio_plot(LikelihoodFit &fit)
//************************************************************************
{
    // --- Getting the current minimizer for likelihood ratio plot ---
    
    TMinuit *minimizer = fit.SetBestFitValues();
   
    // --- Scans parameter 0 ---

    minimizer->Command("SCAn 0");

    // --- Gets the plot from the minimizer and returns it ---

    TGraph *ratio = dynamic_cast<TGraph*>(minimizer->GetPlot());

    
    return ratio;

}

// --- This function sets the best fit values and returns the minimizer ---

//************************************************************************
TMinuit* LikelihoodFit::SetBestFitValues()
//************************************************************************
{
    TMinuit *min = new TMinuit(npar);
    
    best_fit = get_best_fit_values(min);

    return min;

}

// --- This function adds parameters values for the case of fitting with two parameters ---

//************************************************************************
void LikelihoodFit::AddToParameters(double par1, double par2, double par3, double par4, std::string par5)
//************************************************************************
{
    // --- Adding parameters to the vectors ---
    
    par.push_back(par1);
    stepSize.push_back(par2);
    minVal.push_back(par3);
    maxVal.push_back(par4);
    parName.push_back(par5);
}

// --- This function plots the profile likelihood ratio estimation ---

//************************************************************************
TGraph* LikelihoodFit::get_profile_likelihood_ratio_plot(const json &input)
//************************************************************************
{
    // --- Getting the minimum and maximum values for the profile likehood ratio plot ---
    
    double minm = input[region]["mu_min"].get<double>();
	double maxm = input[region]["mu_max"].get<double>();

    // --- 100 points will be used to draw the profile likelihood plot ---
    
    const int npoints = 100;

    // --- 100 values will be used for the mu estimator ---
    // --- So, of course there will be 100 values of the likelihood function ---

    double mu[npoints];
    double lnL[npoints];

    // --- Vectors that will extract the maximum likelihood value for x, for different values of mu ---

    double best_fit_x[npoints];
    double best_fit_x_error[npoints];

    // --- Arrays to be used for the profile likelihood estimation ---

    double numerator[npar];
    double denominator[npar];

    // --- Getting the maximum likelihood values to denominator --

    denominator[0] = best_fit.at(0);
    denominator[1] = best_fit.at(1);

    // --- Defining a TMinuit object ---

    TMinuit *minimizer = new TMinuit(npar);

    // --- Setting the function to be minimized ---
    
    minimizer->SetFCN(LikelihoodFit::fcn);

    for(int i = 0; i < npoints; i++)
    {
        // --- Disecting the maxm - minm interval into 100 points ---
        
        mu[i]= minm+(((maxm-minm)/npoints)*i);

        // --- Initiating the likelihood function value for every loop ---
        
        lnL[i]=0;

        // --- Setting the value of the mu parameter ---

        par[0] = mu[i];

        // --- Giving the parameters to the minimizer ---

        for(int j = 0; j < npar; j++)
        {
            minimizer->DefineParameter(j,parName[j].c_str(),par[j],stepSize[j],minVal[j],maxVal[j]);
        }

        // --- Fixing the mu parameter ---
        
        minimizer->FixParameter(0);

        // --- Minimizing ---

        minimizer->Migrad();

        // --- Getting the maximum likelihood value for x, at a fixed mu ---

        minimizer->GetParameter(1,best_fit_x[i],best_fit_x_error[i]);

        // --- The bin number is the same for every histogram ---

        int nbins = histos.at(0)->GetNbinsX();

        for(int k = 1; k <= nbins; k++)
        {
            numerator[0] = mu[i];
            numerator[1] = best_fit_x[i];

            lnL[i] = lnL[i] + log_likelihood(k, numerator, 1) - log_likelihood(k, denominator, 1);
            
        }

        minimizer->Release(0);

        lnL[i] = -2*lnL[i];
    }

    TGraph *prof_ratio = new TGraph(npoints,mu,lnL);

    return prof_ratio;


}

// --- This function styles the ratios ---

//************************************************************************
void LikelihoodFit::style_ratios(TGraph *plot, int color, const json &input)
//************************************************************************
{
    // --- Getting the X axis limits for the likelihood ratio ---

    double likelihood_ratio_xmin = input[region]["likelihood_ratio_xmin"].get<double>();
    double likelihood_ratio_xmax = input[region]["likelihood_ratio_xmax"].get<double>();
   
    // --- Style the likelihood ratio plot ---

    plot->SetLineWidth(4);
    plot->SetLineColor(color); 
    plot->GetYaxis()->SetRangeUser(0,1.7); 
    plot->GetXaxis()->SetRangeUser(likelihood_ratio_xmin,likelihood_ratio_xmax);
    plot->SetTitle("");
    plot->GetXaxis()->SetTitle("#hat{#mu}");
    plot->GetYaxis()->SetTitle("-2ln#lambda(x;#hat{#mu})");
    plot->SetMinimum(0);
    
}

// --- This function draw the sigma logo on the likelihood plot ---

//************************************************************************
void LikelihoodFit::draw_sigma(double x, double y)
//************************************************************************
{
   TLatex *sigma = new TLatex(x,y,"#color[2]{#pm1#sigma}");
   sigma->SetTextSize(0.05);
   sigma->SetLineColor(2);
   sigma->Draw("same");
}

// --- This function takes the likelihood plots and fits a quadratic function to calculate the errors ---

//************************************************************************
std::tuple<double,double,double> LikelihoodFit::quadratic_fit(TGraph *plot, const json &input)
//************************************************************************
{
    // --- Getting the x values for the fit from the json file ---

    double fit_x1 = input[region]["fit_ratio_x1"].get<double>();
    double fit_x2 = input[region]["fit_ratio_x2"].get<double>();

    // --- Fiting the likelihood function with a second order polynomial around the minimum ---

    TF1 *f1=new TF1("f1","pol2(0)",fit_x1,fit_x2); 
    f1->SetLineColor(kRed); 
    plot->Fit(f1,"R+0Q");

    // --- Retrieving the coefficients of the polynomial ---
    // --- The polynomial is defined as [0] + [1]*x + [2]*x*x

    double a = f1->GetParameter(2);
    double b = f1->GetParameter(1);
    double c = f1->GetParameter(0);

    // --- We want to solve the equation a*x^2 + b*x + c = 1 + minimum ---
    // --- According to the statistical theory, the x values that satisfy this equation correspond to best_fit + 1σ and best_fit - 1σ ---
    // --- This is a typical second order polynomial equation, so we compute the discriminant ---

    double discriminant = pow(b,2)-(4*a*(c-1-f1->GetMinimum()));

    // --- Extracting the solutions of the equation ---

    double x1= (-b+sqrt(discriminant))/(2*a);
    double x2= (-b-sqrt(discriminant))/(2*a);

    double error1 = x1 - f1->GetMinimumX();
    double error2 = f1->GetMinimumX() - x2;

    double total_error = ( error1 + error2 )/2;

    std::cout << f1->GetMinimumX() <<  " +/- " << total_error << std::endl;

    // --- In this point, the user of the code might think that there is a problem with the fit, due to the low value of chi^2/ndf ---
    // --- The answer is that the fitted function (likelihood function), is a smooth function, not a collection of random measurements ---
    // --- Thus, there is no statistical fluctiation associated to each point ---
    // --- In this context, the chi^2/ndf goodness of fit indicator loses its statistical meaning, and a very low value is expected ---
    // --- This is because according to statistical theory, the likelihood function is expected to behave as a parabola close to its minimum ---
    
    std::cout << "Chi squared over NDF: " << f1->GetChisquare()/f1->GetNDF() << std::endl;

    return std::make_tuple(x1,x2,total_error);
}

// --- This function draws the legend for the likelihood comparison plot ---

//************************************************************************
void LikelihoodFit::draw_likelihood_legend(TGraph *ratio, TGraph *profile_ratio)
//************************************************************************
{
   TLegend *leg = new TLegend(0.72,0.15,0.9,0.25);
   leg->AddEntry(ratio,"#bf{Stat}.","l");
   leg->AddEntry(profile_ratio,"#bf{Stat. #oplus Syst.}","l");
   leg->SetBorderSize(0);
   leg->Draw();
}

// --- This function draws the two ratios together ---

//************************************************************************
void LikelihoodFit::draw_ratios(TGraph *ratio, TGraph *profile_ratio,std::unique_ptr<TCanvas> &c, const json &input)
//************************************************************************
{
    // --- Showing the canvas ---

    c->cd();

    // --- Styling the ratios ---

    style_ratios(ratio,kBlue,input);
    
    style_ratios(profile_ratio,kMagenta,input);


    // --- Fiting the two ratios with a quadratic function and printing the information ---

    std::cout << std:: endl << "----------------------------------------------------------------------" << std::endl;

    std::cout << "Input variable: " << var << std::endl;

    std::cout << "----------------------------------------------------------------------" << std::endl;
    
    std::cout << "Likelihood ratio estimation: ";
    
    std::tuple<double,double,double> errors_r = quadratic_fit(ratio,input);

    std::cout << "----------------------------------------------------------------------" << std::endl;
    
    std::cout << "Profile Likelihood ratio estimation: ";
    
    std::tuple<double,double,double> errors_pr = quadratic_fit(profile_ratio,input);

    std::cout << "----------------------------------------------------------------------" << std::endl;

    std::cout << "Luminosity error: " << sqrt(pow(std::get<2>(errors_pr),2) - pow(std::get<2>(errors_r),2)) << std::endl;

    std::cout << "----------------------------------------------------------------------" << std::endl;


    // --- Drawing the likelihood ratio ---

    ratio->Draw("ac");

    // --- Drawing the vertical lines to indicate x1 and x2 for the likelihood ratio ---

    draw_line(kBlue,2,std::get<0>(errors_r),0,std::get<0>(errors_r),1);
    draw_line(kBlue,2,std::get<1>(errors_r),0,std::get<1>(errors_r),1);

    // --- Drawing the profile likelihood ratio ---

    profile_ratio->Draw("c");

    // --- Drawing the vertical lines to indicate x1 and x2 for the profile likelihood ratio ---

    draw_line(kMagenta,2,std::get<0>(errors_pr),0,std::get<0>(errors_pr),1);
    draw_line(kMagenta,2,std::get<1>(errors_pr),0,std::get<1>(errors_pr),1);

    // --- Getting the x value for the sigma logo ---

    double sigma_x = input[region]["sigma_x"].get<double>();

    // --- Drawing the sigma logo ---

    draw_sigma(sigma_x,1.05);

    // --- Getting the x limits for the line ---

    double line_x1 = input[region]["line_xmin"].get<double>();
    double line_x2 = input[region]["line_xmax"].get<double>();

    // --- Draw the common line that intersects the plots at y=1 ---

    draw_line(15,2,line_x1,1,line_x2,1);

    // --- Getting the x value for the atlas logo ---

    double ratio_logo_x = input[region]["ratio_logo_x"].get<double>();

    // --- Draw the ATLAS logo ---

    draw_logo(ratio_logo_x,1.5,1.35,1.2,input);

    // --- Draws the legend ---

    draw_likelihood_legend(ratio,profile_ratio);

}