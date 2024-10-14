# XPRS

eXplainable PRS (XPRS) is a software designed to enhance the interpretability of PRS by decomposing them into gene/region and SNP contribution scores. 
By utilizing Shapley Additive Explanations (SHAP), XPRS calculates the attributed value of each gene or region, providing detailed insights into which genes significantly contribute to an individual’s PRS.

# Getting Started

### 1. Clone the Repository:
XPRS leverages PLINK for certain analyses. The PLINK executable is included in the repository. Before running PLINK, you need to grant execute permissions:
``` 
git clone https://github.com/nayeonkim93/XPRS.git
cd XPRS
 ```
### 3. Plink
XPRS leverages PLINK for certain analyses. The PLINK executable is included in the repository. Before running PLINK, you need to grant execute permissions:
``` 
chmod +x plink
 ``` 
### 3. Install Python and R packages
XPRS is built using both Python and R, and requires specific libraries for proper functionality. To install all required Python and R libraries, run the `install.sh` script. This script will install both the Python dependencies listed in the `requirements.txt` file and the R libraries using the provided `install_packages.R` script.

To execute the `install.sh` script, use the following command:

```
bash install.sh
``` 

### 4. Run XPRS
``` 
export FLASK_APP=app
export FLASK_ENV=development
export FLASK_DEBUG=1
flask run -p 5000
 ```
