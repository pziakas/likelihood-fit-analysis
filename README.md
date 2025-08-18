# Maximum Likelihood Estimation Analysis in Particle Physics Data

## Overview 

This project demonstrates advanced parameter estimation techniques using Maximum Likelihood Fits implemented in C++ with the ROOT Data Analysis Framework. It focuses on estimating key parameters and quantifying their uncertainties, incorporating both statistical and systematic effects through Profile Likelihood Ratio methods. The workflow includes data loading, fitting, and visualization, showcasing proficiency in statistical modeling and data visualization relevant for complex data analysis challenges.

## Features

- **Maximum Likelihood Estimation (MLE):** Implements parameter estimation using ROOT's `TMinuit` minimizer for robust statistical inference.
- **Profile Likelihood Ratio:** Quantifies the impact of nuisance parameters on signal strength estimation, improving model reliability.
- **Uncertainty Quantification:** Fits likelihood curves with quadratic functions to extract parameter uncertainties and confidence intervals.
- **Customizable Analysis Regions:** Supports multiple physics control regions (`ZZ`, `WZ_qcd`), allowing flexible application of the model.
- **Visualization:** Generates publication-quality plots including likelihood ratio curves, profile likelihood comparisons, and annotated confidence intervals.
- **Interactive Workflow:** User inputs guide the analysis, facilitating exploration of different variables and regions.
- **Integration with JSON Configuration:** Uses JSON files to configure input parameters and file paths, enhancing reproducibility and ease of use.

## Technologies Used

- C++17  
- [ROOT](https://root.cern/) (data analysis framework)  
- [JSON for Modern C++](https://github.com/nlohmann/json) (config parsing)

## Usage 

To run the project, the user needs to have the following installed:

- [ROOT](https://root.cern/)
- [JSON for Modern C++](https://github.com/nlohmann/json)

Make sure to configure the path in the JSON file, according you system. Then using the Makefile provided in this project, the user can run the project using the following commands:

```bash 

cd path/to/the/project/
make
./likelihood_fit 
```

The user will need to follow the ineractive prompts to select a variable as input and select a region in which the fit will be performed. 

## Further Reading

For a detailed explanation of the statistical methods and physics context, please refer to my thesis, included in this repository.

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

In case you need further clarification or you are interested in collaborating, feel free to reach out at **panagiotis.ziakas@cern.ch**.
