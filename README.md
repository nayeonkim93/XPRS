# XPRS

eXplainable PRS (XPRS) is a software designed to enhance the interpretability of PRS by decomposing them into gene/region and SNP contribution scores. 
By utilizing Shapley Additive Explanations (SHAP), XPRS calculates the attributed value of each gene or region, providing detailed insights into which genes significantly contribute to an individual’s PRS.

# Dependencies 
To ensure XPRS runs smoothly, several dependencies are required. Below are the steps to set up the necessary software and libraries.
1. Plink
XPRS leverages PLINK for certain analyses. The PLINK executable is included in the repository. Before running PLINK, you need to grant execute permissions:
''' chmod +x plink '''
2. Python libraries
XPRS is built using Python and requires specific libraries to function correctly. All necessary Python dependencies are listed in the requirements.txt file. To install them, execute the following command:
'''pip install -r requirements.txt'''
3. R packages
Some components of XPRS depend on R packages. To install the required R libraries, you can use the provided install_packages.R script. Run the script in your R environment:
'''source("install_packages.R")'''

# Getting Started
