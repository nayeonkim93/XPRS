# Read package names from the text file
packages <- readLines("r_packages.txt")

# Install packages if they are not already installed
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
}
