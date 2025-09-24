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

using json = nlohmann::json;


// --- Declaring histogramms ---

TH1F *h_data;
TH1F *h_ZZ;
TH1F *h_WZ_qcd;
TH1F *h_WZ_ew;
TH1F *h_fakes;
TH1F *h_ttbarV_VVV;
TH1F *h_signal;

// --- Declaring region ---

std::string region;

// --- Declarinng the flag ---

int fit_flag;

// --- Declaring the vectors ----

std::vector <double> best_fit;
std::vector <TH1F*> histos;


// --- Declaring used functions ---


// --- This function plots the line of the ratio plot ---

//************************************************************************
void draw_line(int color, int line_width, double xmin, double ymin, double xmax, double ymax)
//************************************************************************
{
	
	TLine *l = new TLine(xmin,ymin,xmax,ymax);
    l->SetLineColor(color);
	l->SetLineWidth(line_width);
    l->SetLineStyle(9);
    l->Draw("same");
}

// --- This function draws the ATLAS logo ---

//************************************************************************
void draw_logo(double xmin, double xmax1, double xmax2, double xmax3, std::string region, const json &input)
//************************************************************************
{
    TString region_logo(input[region]["region"].get<std::string>());
    
    TLatex *t = new TLatex(xmin,xmax1,"#it{ATLAS} #bf{work in progress}");
    TLatex *t1=new TLatex(xmin,xmax2,"#bf{#sqrt{s}=13 TeV, #intLdt = 139 fb^{-1}}");
    TLatex *t2=new TLatex(xmin,xmax3,region_logo);
    
    t->SetTextSize(0.04);
    t1->SetTextSize(0.035);
    t2->SetTextSize(0.035);
    
    t->Draw();
    t1->Draw();
    t2->Draw();
}


// --- This function draws the legend of the plot ---


//************************************************************************
void draw_legend(const std::vector <TH1F*> &hists, const std::map <TH1F*,std::vector<TString>> &leg_map)
//************************************************************************
{
	TLegend *leg = new TLegend(0.6,0.55,0.85,0.85); 
	for(const auto &hist : hists) leg->AddEntry(hist,leg_map.at(hist).at(0),leg_map.at(hist).at(1));
	leg->SetBorderSize(0);
	leg->Draw();
}

// --- This function ensures that the variable the user provides is valid ---

//************************************************************************
void var_choice(TString &var)
//************************************************************************
{
	std::cout << std::endl << "Enter the variable you want to plot: ";
	std::cin >> var;

	while(var != "mwz" && var != "sum3pt" && var != "mtwz")
	{
		std::cout << "Not valid variable!" << std::endl;
        std::cout << "You can choose between: mwz, sum3pt or mtwz!" << std::endl;
		std::cout << "Enter the variable you want to plot: ";
		std::cin >> var;
	}
}


// --- This function gets the histogram from the file ---

//************************************************************************
void get_hists(const TString &name, TH1F *&h, TFile *file, int fill_color)
//************************************************************************
{
	h = dynamic_cast<TH1F*>(file->Get(name));
 
    if(h == nullptr) std::cout << "Something is wrong" << std::endl;

	style_hists(h,fill_color);
}

// --- This function styles the histograms ---

//************************************************************************
void style_hists(TH1F *&h, int fill_color)
//************************************************************************
{
    TString name = h->GetName();
	
    if(name.Contains("data"))
	{
		h->SetMarkerStyle(20);
		h->SetMarkerColor(kBlack);
		h->SetLineColor(kBlack);
	}

	else h->SetFillColor(fill_color);
}

//************************************************************************
void draw_ratio_plot(const TString &var, TCanvas *&c, const std::vector <TH1F*> &hists, TH1F *&h_data,const std::vector <TH1F*> &hists_leg, const std::map <TH1F*,std::vector<TString>> &leg_map,std::string region, const json &input)
//************************************************************************
{
	int size = hists.size();

    // --- Converting the TStirng var into a string object ---

    std::string string_var(var.Data());

    // --- Defining the stack ---

     THStack *s = new THStack("stack",";; Counts");
	
	// --- Getting the stack limits from the json file ---

	double stack_max = input[region][string_var]["stack_upper_lim"].get<double>();
    double stack_min = input[region][string_var]["stack_lower_lim"].get<double>();


	// --- Filling the stack ---

	for(const auto &hist : hists) s->Add(hist);

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

    double logo_y1 = input[region][string_var]["atlas_logo_y1"].get<double>();
    double logo_y2 = input[region][string_var]["atlas_logo_y2"].get<double>();
    double logo_y3 = input[region][string_var]["atlas_logo_y3"].get<double>();
    double logo_x = input[region]["stack_logo_x"].get<double>();

    draw_logo(logo_x,logo_y1,logo_y2,logo_y3,region,input);

    // --- Draw the legend ---
    
    draw_legend(hists_leg,leg_map);

    // --- Draw data ---

    h_data->Draw("e1x0p same");

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

    TH1F *h_Data = dynamic_cast<TH1F*>(h_data->Clone("h_Data"));
    TH1F *h_sum = dynamic_cast<TH1F*>(hists.at(0)->Clone("h_sum"));

    // --- Not taking again into account the first histo ---
    // --- This is because it was already copied above ---
    // --- Otherwise its contribution is counted twice ---

    for(int i = 1; i < size; i++) h_sum->Add(hists.at(i));

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

    int first_bin = h_data->FindFirstBinAbove(-10);
    int last_bin = h_data->FindLastBinAbove(-10);

    double minx = h_data->GetBinLowEdge(first_bin);
    double maxx = h_data->GetBinLowEdge(last_bin) + h_data->GetBinWidth(last_bin);

    draw_line(kBlack,1,minx,1,maxx,1);
 }

 // --- This function gives the binned likelihood value ---
 // --- The value that is returned is different for each bin ---
 // --- The flag corresponds to the calculation with or withour the nuisance parameter ---

 //************************************************************************
double log_likelihood(int bin, double par[], int flag)
//************************************************************************
{
    double val;

    TF1 *poisson = new TF1("p","TMath::Poisson(x,[0])",0,1e5);
    TF1 *gauss = new TF1("g","gaus(0)",-10,10);

    // --- Getting the singal counts from the global signal histo ---

    double s = h_signal->GetBinContent(bin);

    // --- Getting the total background counts ---

    double b = 0;

    for(auto &hist : histos) b += hist->GetBinContent(bin);

    // --- The luminosity systematic is inserted in the likelihood as a nuisance parameter p[0] which can change the number of total predicted events ---
    // --- This parameter is known to follow a Gaussian distribution, is centered around 1 and has a realative error of 2% ---

    gauss->SetParameter(0,1/(sqrt(TMath::TwoPi())*(0.02)));
    gauss->SetParameter(1,1);
    gauss->SetParameter(2,0.02);

    // --- Corresponds to the non-systematics case log-likelihood estimation ---
 
    if(flag == 0)
    {
        poisson->SetParameter(0,( (par[0]*s) + b ) );
       
        val = log(poisson->Eval(h_data->GetBinContent(bin)));
    } 

    // --- Corresponds to the calculation with the systematic uncertainty ---

    else if(flag == 1)
    {
        poisson->SetParameter(0,par[1]*( (par[0]*s) + b ) );

        val = log(poisson->Eval(h_data->GetBinContent(bin))) + log(gauss->Eval(par[1]));
    } 

    return val;

}


// --- This is the function that is going to be given to the minimizer ---
// --- It can't get different arguments, otherwise TMinuit is not working ---

//************************************************************************
void fcn(int &npar, double *deriv, double &f, double *par, int flag)
//************************************************************************
{
    double lnL = 0;

    int size = best_fit.size();

    double vals[size];

    for(int i = 0; i < size; i++) vals[i] = best_fit.at(i);

    // --- All of the histograms have the same number of bins ---

    int nbins = h_data->GetNbinsX();

    for(int i = 1; i <= nbins; i++)
    {
        // --- Binned likelihood fit ---

        if(fit_flag == 0)
        {
            lnL = lnL + log_likelihood(i, par, 0);
        }

        // --- Binned likelihood ratio fit ---

        else if(fit_flag == 1) 
        {
            vals[0] = best_fit.at(0);
            
            lnL = lnL + log_likelihood(i, par, 0) - log_likelihood(i, vals, 0);
        }

         // --- Likelihood ratio fit with one nuisance parameter ---

        else if(fit_flag == 2)
        {
            lnL = lnL + log_likelihood(i, par, 1);
        } 
    }

    f = -2*lnL;
}

// --- This function gets a TMinuit minimizer of npar parameters and returns a vector of the best fit values ---

//************************************************************************
std::vector<double> get_best_fit_values(TMinuit *&minimizer,const int npar, const std::vector <double> &par, const std::vector <double> &stepSize, const std::vector <double> &minVal,const std::vector <double> &maxVal, const std::vector <std::string> &parName)
//************************************************************************
{
    std::vector<double> values;

    minimizer->SetFCN(fcn);

    for(int i=0;i<npar;i++) minimizer->DefineParameter(i,parName[i].c_str(),par[i],stepSize[i],minVal[i],maxVal[i]);

    
    // --- The algorithm that is called for the minimization ---
    
    minimizer->Migrad();
   
    double outpar[npar];
    double error[npar];

   for(int i=0;i<npar;i++)
   {
      minimizer->GetParameter(i,outpar[i],error[i]);

      std::cout << parName[i] << ": " << outpar[i] << " +/- " << error[i] << std::endl;

      values.push_back(outpar[i]);
    
    }


   return values;

   
}

// --- This function plots the likelihood ratio estimation ---

//************************************************************************
TGraph* get_likelihood_ratio_plot(TMinuit *&minimizer)
//************************************************************************
{

    // --- Scans parameter 0 ---

    minimizer->Command("SCAn 0");

    
    // --- Gets the plot from the minimizer and returns it ---

    TGraph *ratio = dynamic_cast<TGraph*>(minimizer->GetPlot());

    return ratio;
    
}

// --- This function draw the sigma logo on the likelihood plot ---

//************************************************************************
void draw_sigma(double x, double y)
//************************************************************************
{
   TLatex *sigma = new TLatex(x,y,"#color[2]{#pm1#sigma}");
   sigma->SetTextSize(0.05);
   sigma->SetLineColor(2);
   sigma->Draw("same");
}

// --- This function styles the ratios ---

//************************************************************************
void style_ratios(TGraph *plot, int color, std::string region, const json &input)
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

// --- This function plots the profile likelihood ratio estimation ---

//************************************************************************
TGraph* get_profile_likelihood_ratio_plot(TMinuit *&minimizer, double minm, double maxm, const int npars, std::vector <double> &par, const std::vector <double> &stepSize, const std::vector <double> &minVal,const std::vector <double> &maxVal, const std::vector <std::string> &parName)
//************************************************************************
{
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

    double numerator[npars];
    double denominator[npars];

    // --- Getting the maximum likelihood values to denominator --

    denominator[0] = best_fit.at(0);
    denominator[1] = best_fit.at(1);

    for(int i = 0; i < npoints; i++)
    {
        // --- Disecting the maxm - minm interval into 100 points ---
        
        mu[i]= minm+(((maxm-minm)/npoints)*i);

        // --- Initiating the likelihood function value for every loop ---
        
        lnL[i]=0;

        // --- Setting the value of the mu parameter ---

        par[0] = mu[i];

        // --- Giving the parameters to the minimizer ---

        for(int j = 0; j < npars; j++)
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

        int nbins = h_data->GetNbinsX();

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

// --- This function takes the likelihood plots and fits a quadratic function to calculate the errors ---

//************************************************************************
std::tuple<double,double,double> quadratic_fit(TGraph *plot, std::string region, const json &input)
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

// --- This function draws the two ratios together ---

//************************************************************************
void draw_ratios(TGraph *ratio, TGraph *profile_ratio, TCanvas *canvas, const TString &var, std::string region, const json &input)
//************************************************************************
{
    // --- Showing the canvas ---

    canvas->cd();

    // --- Styling the ratios ---

    style_ratios(ratio,kBlue,region,input);
    
    style_ratios(profile_ratio,kMagenta,region,input);


    // --- Fiting the two ratios with a quadratic function and printing the information ---

    std::cout << std:: endl << "----------------------------------------------------------------------" << std::endl;

    std::cout << "Input variable: " << var << std::endl;

    std::cout << "----------------------------------------------------------------------" << std::endl;
    
    std::cout << "Likelihood ratio estimation: ";
    
    std::tuple<double,double,double> errors_r = quadratic_fit(ratio,region,input);

    std::cout << "----------------------------------------------------------------------" << std::endl;
    
    std::cout << "Profile Likelihood ratio estimation: ";
    
    std::tuple<double,double,double> errors_pr = quadratic_fit(profile_ratio,region,input);

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

    draw_logo(ratio_logo_x,1.5,1.35,1.2,region,input);

    // --- Draws the legend ---

    draw_likelihood_legend(ratio,profile_ratio);

}

// --- This function draws the legend for the likelihood comparison plot ---

//************************************************************************
void draw_likelihood_legend(TGraph *ratio, TGraph *profile_ratio)
//************************************************************************
{
   TLegend *leg = new TLegend(0.72,0.15,0.9,0.25);
   leg->AddEntry(ratio,"#bf{Stat}.","l");
   leg->AddEntry(profile_ratio,"#bf{Stat. #oplus Syst.}","l");
   leg->SetBorderSize(0);
   leg->Draw();
}

// --- This function is asking the user if they want to print the canvas in a pdf form ---

//************************************************************************
void print_canvas(TCanvas *&c, TString image_name)
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
        c->Print(image_name + ".pdf");
    }

}

// --- This function ensures that the region the user provides is valid ---

//************************************************************************
void region_choice(std::string &region)
//************************************************************************
{
    std::cout << std::endl << "Enter the region you want to plot: ";
    std::cin >> region;

    while(region != "ZZ" && region != "WZ_qcd")
    {
        std::cout << "Not valid region!" << std::endl;
        std::cout << "You can choose between: ZZ and WZ_qcd!" << std::endl;
        std::cout << "Enter the region you want to plot: ";
        std::cin >> region;
    }
}

// --- This function finds the histogram that corresponds to the signal in each region ---
// --- After that it removes it from the global vector, so that all the other histograms correspond to backgrounds ---

//************************************************************************
void find_signal_hist(std::vector<TH1F*> &histos)
//************************************************************************
{
    // --- Initializing the for loop using an iterator to the vector ---

    for(auto it = histos.begin(); it != histos.end(); )
    {
        TString name = (*it)->GetName();
        bool condition = name.Contains(region);

        // --- If the condition is satisfied, the loop ends ---
        // --- The signal histogram is removed from the vector ---
        
        if(condition)
        {
            h_signal = (*it);
            it = histos.erase(it);
            break;
        }

        else it++;
    }
}


