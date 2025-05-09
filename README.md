# XPRS
eXplainable PRS (XPRS) is a software designed to enhance the interpretability of PRS by decomposing them into gene/region and SNP contribution scores. 
By utilizing Shapley Additive Explanations (SHAP), XPRS calculates the attributed value of each gene or region, providing detailed insights into which genes significantly contribute to an individual’s PRS. 
For more information, see our manuscript 

[XPRS: A Tool for Interpretable and Explainable Polygenic Risk Score](https://doi.org/10.1093/bioinformatics/btaf143)

## Access the XPRS Web Service

You can explore XPRS via our cloud-based web service:  
[https://xprs.leelabsg.org](https://xprs.leelabsg.org)

### How to View Preprocessed Data:

To view the preprocessed data, follow these steps:

1. **Navigate to the `Run` section** and select **Case 1-1**.
2. **Enter the following paths** in the respective fields:
   - **Cohort.genotype.file**: `./data/sample`
   - **RDS.file**: `./output/sample.RDS`
3. **Submit the form** to view the results.

For **data security reasons**, you cannot upload data directly through the web interface. Therefore you must run XPRS locally. 

## Running XPRS Locally

To use XPRS, you must run it locally and specify the direct paths to your data files on your machine. 

To run XPRS locally:

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
chmod +x start_rserve.sh
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

### 5. When running XPRS, make sure to specify the **local paths** to your files. 

For example:

```
file: /path/to/your/data/sample
```
