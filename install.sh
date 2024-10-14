# install.sh
#!/bin/bash

# Install R packages
Rscript install_packages.R

# Install Python packages
pip install -r requirements_py.txt
