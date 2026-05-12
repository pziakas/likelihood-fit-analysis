#include "likelihood_fit.h"

int main(int argc, char **argv)
{
	auto c = std::make_unique<TCanvas>("c", "Canvas", 800, 800);

    LikelihoodFit fit;

    json input = fit.LoadJSON();

    try
    {
        fit.load_data(fit,input);

        fit.visualize_data(fit,"distribution",c,input);

        fit.perform_fit(fit,c,input);
    }   

    catch(const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}