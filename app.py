import os
import sqlite3
import subprocess
import glob
import atexit
from flask import Flask, render_template, request, url_for, flash, redirect, Response
from werkzeug.exceptions import abort
import pyRserve
import gzip
import shutil

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your secret key'

# Path to the gzipped file
gzipped_file_path = './data/cS2G_annotation.txt.gz'
output_file_path = './data/cS2G_annotation.txt'

# Check if the gzipped file exists and the unzipped file does not
if os.path.exists(gzipped_file_path) and not os.path.exists(output_file_path):
    with gzip.open(gzipped_file_path, 'rb') as f_in:
        with open(output_file_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f'Unzipped {gzipped_file_path} to {output_file_path}')


# Function to check if Rserve is running
def is_rserve_running():
    try:
        result = subprocess.run(['pgrep', '-f', 'Rserve'], stdout=subprocess.PIPE)
        return result.returncode == 0
    except subprocess.CalledProcessError:
        return False


# Start Rserve if it's not already running
def start_rserve():
    try:
        # Check if Rserve is already running
        result = subprocess.run(['pgrep', '-f', 'Rserve'], stdout=subprocess.PIPE)
        if result.returncode != 0:  # Rserve is not running
            subprocess.run(['./start_rserve.sh'], check=True)
        else:
            print("Rserve is already running.")
    except subprocess.CalledProcessError as e:
        print("Error starting Rserve: ", e)

# Ensure Rserve is running
start_rserve()

# Connect to the Rserve instance
try:
    conn = pyRserve.connect()
except pyRserve.rexceptions.RConnectionRefused as e:
    print("Could not connect to Rserve:", e)
    exit(1)

# Determine the absolute path to the project directory
project_dir = os.path.abspath(os.path.dirname(__file__))
print(f"Project directory: {project_dir}")

# Replace backslashes with forward slashes for R compatibility (if needed)
project_dir = project_dir.replace("\\", "/")
print(f"Formatted project directory for R: {project_dir}")


# Set the working directory and load necessary libraries and data once
conn.eval(f'''
    setwd("{project_dir}")
    library(parallel)
    library(dplyr)
    library(readr)
    library(data.table)
    library(hash)
    library(optparse)
    library(Rcpp)
    library(parallel)
    library(ggplot2)
    library(ggrepel)
    library(colorspace)
    library(plotly)
    library(htmltools)
    library(htmlwidgets)
    library(jsonlite)
    source("R/which_symbol_important_individual.R")
    source("R/which_symbol_important_population.R")
    source("R/h2_cut.R")
    source("R/map_snps.R")
    source("R/map_snps_no_gwas.R")
    sourceCpp("./C++/BedFileReader.cpp")

    setClass("BedFileReader", representation( pointer = "externalptr" ) )
    BedFileReader_method <- function(name) {{paste( "BedFileReader", name, sep = "__" ) }}
    setMethod( "$", "BedFileReader", function(x, name ) {{function(...) .Call(BedFileReader_method(name),x@pointer, ... )}})
    setMethod("initialize","BedFileReader", function(.Object, ...) {{
    .Object@pointer <-.Call(BedFileReader_method("new"), ... )
    .Object}})
''')

# Define the function to load data if it exists
conn.eval('''
load_data_if_exists <- function(file_path) {
  cat("Attempting to load RDS file from:", file_path, "\n")
  if (file.exists(file_path)) {
    data <- readRDS(file_path)
    # Check for the presence of the components
    if (all(c("annotated_info", "snp_info", "important_genes", "full_prs", "window_size") %in% names(data))) {
      annotated_info <<- data$annotated_info
      snp <<- data$snp_info
      important_genes <<- data$important_genes
      full_PRS <<- data$full_prs
      write.table(data$window_size, file = "./static/data/window_size.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

      cat("Data loaded successfully with annotated_info and snp_info, important genes, full_PRS and window_size from", file_path, "\n")
      
      # Check if full_PRS and PRS_info are also present
      if (all(c("PRS_info") %in% names(data))) {
        PRS_info <<- data$PRS_info
        cat("Data loaded successfully with PRS_info from", file_path, "\n")
       }
     }
    else {
      cat("Required data annotated_info, snp_info, important_genes, full_prs, window_size not found in", file_path, "\n")
    }
   } 
   else {
    cat("File", file_path, "does not exist. No data loaded.\n")
  }
}
''')

# Call the function to load the data
conn.eval('load_data_if_exists("./output/data.RDS")')

# Define a function to read test file if exists, else use default
def get_test_file_path():
    test_file_path = './static/data/test_file.txt'
    if os.path.exists(test_file_path):
        with open(test_file_path, 'r') as file:
            test = file.read().strip()
    else:
        test = './data/new_test'
    return test


def delete_all_files():
    # Define folders and file patterns to delete
    paths_to_delete = [
        './static/data/*',
        './static/data/debug/*',
        './output/data.RDS',
        './data/new_test.*',
        './data/final_indi.*',
        './output/indi_data.RDS',
        'reading_file.log'
    ]
    
    for path in paths_to_delete:
        for file_path in glob.glob(path):
            if os.path.isfile(file_path):
                os.remove(file_path)
                print(f"Deleted {file_path}")

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/tutorial', methods=('GET', 'POST'))
def tutorial():
    return render_template('tutorial.html')

@app.route('/tutorial_input', methods=('GET', 'POST'))
def tutorial_input():
    return render_template('tutorial_input.html')

@app.route('/tutorial_process', methods=('GET', 'POST'))
def tutorial_process():
    return render_template('tutorial_process.html')

@app.route('/run', methods=('GET', 'POST'))
def run():
    return render_template('run.html')

@app.route('/link', methods=('GET', 'POST'))
def link():
    return render_template('link.html')


@app.route('/run_1', methods=('GET', 'POST'))
def run_1():
    if request.method == 'POST':
        test = request.form['test']
        GWAS = request.form['GWAS']
        beta = request.form['beta']
        cS2G = request.form['cS2G']
         # annot = request.form['annot']
        important_genes = request.form['important_genes']

        window = request.form.get('window', 200000)  # Default value if not provided
        with open('./static/data/window_size.txt', 'w') as f:
            f.write(window)
        cpu = request.form['cpu']
        cut = request.form['cut']
        

        if not test:
            flash('population.genomic.file is required!')
    
        if not beta:
            flash('PRS.scoring.file is required!')
                
        if not important_genes:
            flash('important_genes.file is required!')
        
        if test and beta: 
            # with open('./output/test_file.txt', 'w') as f:
            #     f.write(test)
            if not GWAS: 
                if not cut:  
                    if not cS2G:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --number.cpu " + cpu)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window )     
                        else: 
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window  + " --number.cpu " + cpu)     
                    else:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --number.cpu " + cpu)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window )     
                        else:   
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window  + " --number.cpu " + cpu)     
                    
                    return redirect(url_for("result"))
                if cut: 
                    if not cS2G:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --h2.cut " + cut)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --number.cpu " + cpu + " --h2.cut " + cut)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window + " --h2.cut " + cut)     
                        else: 
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window  + " --number.cpu " + cpu + " --h2.cut " + cut)     
                    else:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --h2.cut " + cut)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --number.cpu " + cpu + " --h2.cut " + cut)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window + " --h2.cut " + cut)     
                        else:   
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window  + " --number.cpu " + cpu + " --h2.cut " + cut)     
                    
                    return redirect(url_for("result"))

            elif GWAS: 
                if not cut:  
                    if not cS2G:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --number.cpu " + cpu)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window )     
                        else: 
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window  + " --number.cpu " + cpu)     
                    else:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --number.cpu " + cpu)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window )     
                        else:   
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window  + " --number.cpu " + cpu)     
                    return redirect(url_for("result"))
                if cut: 
                    if not cS2G:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --h2.cut " + cut)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --number.cpu " + cpu + " --h2.cut " + cut)  
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window + " --h2.cut " + cut)     
                        else: 
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window  + " --number.cpu " + cpu + " --h2.cut " + cut)    
                    else:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --h2.cut " + cut)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --number.cpu " + cpu + " --h2.cut " + cut)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window + " --h2.cut " + cut)     
                        else:   
                            os.system("Rscript pre_main_processing.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window  + " --number.cpu " + cpu + " --h2.cut " + cut)     
                    return redirect(url_for("result"))
     
        return render_template('run_1.html')
    return render_template('run_1.html')

@app.route('/run_already_done', methods=('GET', 'POST'))
def run_already_done():
    if request.method == 'POST':
        test = request.form['test']
        rds = request.form['rds']

        if not test:
            flash('population.genomic.file is required!')
    
        if not rds:
            flash('rds.file is required!')
                
        if test and rds:

            test_file_path = './static/data/test_file.txt'
            with open(test_file_path, 'w') as f:
                f.write(test)

            # Pass the test and rds file paths to R and run the processing script
            try:
                conn.r.test_file = test
                conn.r.rds_file = rds
                # Execute the R script and capture the time taken
                time_taken = conn.eval(f'''
                    sink(file = './static/data/debug/run_already_done.txt')
                    start_time <- Sys.time()
                    print(test_file)
                    print(rds_file)

                    # Load the data from the RDS file
                    data <- readRDS(rds_file)
                    annotated_info <- data$annotated_info
                    important_genes <- data$important_genes
                    snp <- data$snp
                    full_PRS <<- data$full_prs
                    PRS_info <<- data$PRS_info
                    write.table(data$window_size, file = "./static/data/window_size.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

                    # To load the image from the RDS file
                    loaded_image_binary <- data$table
                    writeBin(loaded_image_binary, "./static/data/population_gene.png")
                    # To load the html from the RDS file
                    loaded_html <- data$plot
                    htmltools::save_html(html = loaded_html, file = "./static/data/population_gene.html")

                    end_time <- Sys.time()
                    time_taken <- end_time - start_time
                    print(paste0("Time it takes for loading cohort rds_file ",time_taken, " seconds"))
                    sink()
                ''')
            except pyRserve.rexceptions.REvalError as e:
                flash(f"An error occurred while running the R script: {e}")

            return redirect(url_for("result"))

    return render_template('run_already_done.html')

@app.route('/run_1/result', methods=('GET', 'POST'))
def result():
    if request.method == 'POST':
        test = get_test_file_path()
        iid = request.form['iid']

        if not test:
            flash('test.file is required!')
        if not iid:
            flash('IID is required!')
        
        else:
            # Pass the arguments to R and run the analysis
            conn.r.test_file = test
            conn.r.iid = iid
   
            time_taken = conn.eval(f'''
                sink(file = './static/data/debug/run1_result.txt')
    

                # Check if any of the key data objects do not exist or are NULL
                if (!exists("annotated_info", envir = .GlobalEnv) || is.null(annotated_info) ||
                    !exists("important_genes", envir = .GlobalEnv) || is.null(important_genes) ||
                    !exists("snp", envir = .GlobalEnv) || is.null(snp) ||
                    !exists("full_PRS", envir = .GlobalEnv) || is.null(full_PRS) ||
                    !exists("PRS_info", envir = .GlobalEnv) || is.null(PRS_info)) {{

                    print("Here - loading data")
                    load_data_if_exists("./output/data.RDS")
                }}


                start_time <- Sys.time()
        
                print(test_file)
                print(iid)

                which_symbol_is_important_individual(iid, test_file, snp, full_PRS, PRS_info, annotated_info, important_genes,'./static/data/snp_level_manhattan_plot.json')
                end_time <- Sys.time()
                print(paste0("Time it takes for which_symbol_is_important_individual function is ", end_time - start_time, " seconds"))
                sink()
            ''')
       
            
            
            return redirect(url_for("result_individual"))
                
    return render_template('result.html')


@app.route('/run_2', methods=('GET', 'POST'))
def run_2():
    if request.method == 'POST':
        test = request.form['test']
        GWAS = request.form['GWAS']
        beta = request.form['beta']
        cS2G = request.form['cS2G']
        important_genes = request.form['important_genes']

        window = request.form.get('window', 100000)  # Default value if not provided
        with open('./static/data/window_size.txt', 'w') as f:
            f.write(window)
        cpu = request.form['cpu']
        cut = request.form['cut']
        
        if not test:
            flash('reference.genome.file is required!')
        
        if not beta:
            flash('PRS.scoring.file is required!')

        if not important_genes:
            flash('important_genes.file is required!')
        
        if test and beta: 
      
            if not GWAS: 
                if not cut:  
                    if not cS2G:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --number.cpu " + cpu)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window )     
                        else: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window  + " --number.cpu " + cpu)     
                    else:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --beta.file " + beta +" --important_genes.file " + important_genes +  " --cS2G.file " + cS2G + " --number.cpu " + cpu)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window )     
                        else:   
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window  + " --number.cpu " + cpu)     
                    
                    return redirect(url_for("result_2"))
                if cut: 
                    if not cS2G:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --h2.cut " + cut)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --number.cpu " + cpu + " --h2.cut " + cut)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window + " --h2.cut " + cut)     
                        else: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window  + " --number.cpu " + cpu + " --h2.cut " + cut)     
                    else:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --h2.cut " + cut)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --number.cpu " + cpu + " --h2.cut " + cut)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing_2.R --test.file " + test +  " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window + " --h2.cut " + cut)     
                        else:   
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window  + " --number.cpu " + cpu + " --h2.cut " + cut)     
                    
                    return redirect(url_for("result_2"))

            elif GWAS: 
                if not cut:  
                    if not cS2G:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --number.cpu " + cpu)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window )     
                        else: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window  + " --number.cpu " + cpu)     
                    else:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --number.cpu " + cpu)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window )     
                        else:   
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window  + " --number.cpu " + cpu)     
                    return redirect(url_for("result_2"))
                if cut: 
                    if not cS2G:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --h2.cut " + cut)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --number.cpu " + cpu + " --h2.cut " + cut)  
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window + " --h2.cut " + cut)     
                        else: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --window.size " + window  + " --number.cpu " + cpu + " --h2.cut " + cut)    
                    else:
                        if not window and not cpu:
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --h2.cut " + cut)
                        elif not window and cpu: 
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --number.cpu " + cpu + " --h2.cut " + cut)    
                        elif window and not cpu:        
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window + " --h2.cut " + cut)     
                        else:   
                            os.system("Rscript pre_main_processing_2.R --test.file " + test + " --GWAS.file " + GWAS + " --beta.file " + beta + " --important_genes.file " + important_genes + " --cS2G.file " + cS2G + " --window.size " + window  + " --number.cpu " + cpu + " --h2.cut " + cut)     
                    return redirect(url_for("result_2"))
        return render_template('run_2.html')
    return render_template('run_2.html')

@app.route('/run_already_done_2', methods=('GET', 'POST'))
def run_already_done_2():
    if request.method == 'POST':
        test = request.form['test']
        rds = request.form['rds']

        if not test:
            flash('reference.genomic.file is required!')
    
        if not rds:
            flash('rds.file is required!')
                
        if test and rds:
            # Pass the test and rds file paths to R and run the processing script

      
            try:
                conn.r.test_file = test
                conn.r.rds_file = rds
                # Execute the R script and capture the time taken
                time_taken = conn.eval(f'''
                    start_time <- Sys.time()

                    # Load the data from the RDS file
                    sink(file = './static/data/debug/run_already_done_2.txt')
                    print(test_file)
                    print(rds_file)
                
                    data <- readRDS(rds_file)
                    annotated_info <- data$annotated_info
                    snp <- data$snp_info
                    important_genes <- data$important_genes
                    full_PRS <- data$full_prs
                    write.table(data$window_size, file = "./static/data/window_size.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

                    # To load the image from the RDS file
                    loaded_image_binary <- data$table
                    writeBin(loaded_image_binary, "./static/data/population_gene.png")
                    loaded_html <- data$plot
                    htmltools::save_html(html = loaded_html, file = "./static/data/population_gene.html")

                    end_time <- Sys.time()
                    time_taken <- end_time - start_time
                    print(paste0("Time it takes for loading reference rds_file ",time_taken, " seconds"))
                    sink()
                ''')
    
            except pyRserve.rexceptions.REvalError as e:
                flash(f"An error occurred while running the R script: {e}")

            return redirect(url_for("result_2"))

    return render_template('run_already_done_2.html')



@app.route('/run_2/result_2', methods=('GET', 'POST'))
def result_2():
    if request.method == 'POST':
        indi = request.form['indi']
        iid = request.form['iid']
        # window = request.form.get('window', 100000)  # Default value if not provided
        cpu = request.form.get('cpu', 8)  # Default value if not provided

        if not indi:
            flash('individual.genomic.file is required!')
        
        if not iid:
            flash('IID is required!')
        else:
           
            # Pass the arguments to R and run the analysis
            conn.r.indi_file = indi
            conn.r.iid = iid

            # Ensure window and cpu are properly set with defaults if empty

            if not cpu:
                cpu = 8

            conn.r.num_cores = cpu

    

            conn.eval(f'''
                # Redirect output to a file
                sink(file = './static/data/debug/run_2_result_2.txt')


                # Check if any of the key data objects do not exist or are NULL
                if (!exists("annotated_info", envir = .GlobalEnv) || is.null(annotated_info) ||
                    !exists("important_genes", envir = .GlobalEnv) || is.null(important_genes) ||
                    !exists("snp", envir = .GlobalEnv) || is.null(snp) ||
                    !exists("full_PRS", envir = .GlobalEnv) || is.null(full_PRS) ) {{

                    print("Here - loading data")
                    load_data_if_exists("./output/data.RDS")
                }}
            
                start_time_T <- Sys.time()
                write.table(snp$snpid, file = "./data/snp_indi.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
                cmd <- paste0("./plink --bfile ", indi_file, " --extract ./data/snp_indi.txt --make-bed --out ./data/temp_indi")
                system(cmd)

                system("cut -f2 ./data/temp_indi.bim | sort | uniq -d > ./data/duplicates.txt")
                cmd <- paste0("./plink --bfile ./data/temp_indi --exclude ./data/duplicates.txt --make-bed --out ./data/temp_indi_no_dup")
                system(cmd)

                bim <- fread("./data/temp_indi_no_dup.bim")


                # merge the bim file and beta file then match the allele A1/A2
                merge <- left_join(bim, snp,  by = c("V2" = "snpid"))
                merge <- merge[is.na(merge$chr) == F,]
                merge <- merge[(merge$V1 == merge$chr),]
                merge <- merge[(merge$V4 == merge$pos),]

                # if bim A1 == PRS beta A2, multiply -1 
                merge$flip <- ifelse(merge$V5 == merge$A2 & merge$V6 == merge$A1, 1, 
                                        ifelse(merge$V5 == merge$A1 & merge$V6 == merge$A2, 0, NA))

                merge <- merge[!is.na(merge$flip), ]

                # save beta file and return it 
                merge <- merge[,c("V1", "V2", "V4", "A1", "A2", "beta", "original_beta", "AF", "flip", "h2")]
                names(merge) <- c("chr", "snpid", "pos", "A1", "A2", "beta", "original_beta","AF","flip", "h2")


                write.table(merge$snpid, file = "./data/snp_cut_off.txt", row.names = F, col.names=F, quote = F)
                cmd = paste0("./plink --bfile ./data/temp_indi_no_dup --extract ./data/snp_cut_off.txt --make-bed --out ./data/final_indi")
                system(cmd)

                snp <- merge
                indi_file <- "./data/final_indi"

                system("rm -rf ./data/snp_indi.txt")
                system("rm -rf ./data/snp_cut_off.txt")
                system("rm -rf ./data/temp_indi.*")
                system("rm -rf ./data/temp_indi_no_dup.*")
                system("rm -rf ./data/duplicates.txt")

                start_time <- Sys.time()
                
                print(indi_file)
                fam <- fread(paste0(indi_file, ".fam"))
                bim <- fread(paste0(indi_file, ".bim"))
                BedFileReader <- new("BedFileReader", paste0(indi_file, ".fam"), paste0(indi_file, ".bim"), paste0(indi_file, ".bed"))
                indi_PRS <- BedFileReader$calculateWholePRS(snp$snpid, snp$original_beta, snp$flip)
                indi_PRS <- indi_PRS[which(fam$V2 == iid)] 
                end_time <- Sys.time()
                print(paste0("Time it takes to calculate PRS for individual is ", end_time - start_time, " seconds"))
                
                INDI_PRS <- data.frame(IID = iid, PRS = indi_PRS) 

                if (!is.data.frame(full_PRS)) {{
                    full_PRS <- data.frame(IID = rep(NA, length(full_PRS)), PRS = full_PRS)
                }}

                full_PRS <- rbind(full_PRS, INDI_PRS)
                # snp$beta <- snp$original_beta * sd(full_PRS$PRS)
                full_PRS$PRS <- (full_PRS$PRS - mean(full_PRS$PRS)) / sd(full_PRS$PRS) 
                full_PRS <- as.data.frame(full_PRS)
                write.table(full_PRS[full_PRS$IID == iid,]$PRS, file = paste0('./static/data/prs.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
               
                
                snp_map <- data.table(
                snpid = snp$snpid,
                beta = snp$beta,
                flip = snp$flip,  # Ensure this column exists
                AF = snp$AF
                )
                setkey(snp_map, snpid)

                # Ensure the sequence is an integer sequence
                num_rows <- nrow(annotated_info)
                sequence <- seq(1, num_rows)

                start_time <- Sys.time()
                PRS_info <- mclapply(sequence, function(i) {{
                    snp_ids <- annotated_info$snps_included[[i]]
        
                    # Retrieve beta and flip values for these SNPs from snp_map
                    snp_subset <- snp_map[J(snp_ids), .(beta, flip, AF), nomatch = 0]
    

                    snp_ids <- snp_ids[snp_ids %in% snp_map$snpid]

                     # Check if the number of snp_ids matches the number of rows in snp_subset
                    if (length(snp_ids) != nrow(snp_subset)) {{
                        # Log a warning if there is a mismatch
                        print(paste("Mismatch in SNP count for row", i, ": Expected", length(snp_ids), "but found", nrow(snp_subset)))
                    }}

                    # Handle cases where no SNPs are found
                    if (nrow(snp_subset) == 0) {{
                        return(0)  # Or another appropriate default value
                    }}

                    # Extract beta and flip lists
                    betas <- snp_subset$beta
                    flips <- snp_subset$flip
                    AFs <- snp_subset$AF


                    # Call the C++ function to calculate PRS with flipping
                    PRS <- BedFileReader$calculatePRS(snp_ids, betas, flips, AFs)

                    
                    return(PRS)
                    }}, mc.cores = num_cores, mc.set.seed = TRUE)
                    end_time <- Sys.time()
                    print(paste0("Time it takes to calculate Gene Contribution Score for individual is ", end_time - start_time, " seconds"))


                start_time <- Sys.time()
                PRS_info_matrix <- do.call(cbind, PRS_info)
                # PRS_info_matrix[which(fam$V2 == iid), ] <- PRS_info_matrix[which(fam$V2 == iid), ] - (annotated_info$mean * (nrow(full_PRS) -1) + PRS_info_matrix[which(fam$V2 == iid), ])/nrow(full_PRS)
                PRS_info_matrix[which(fam$V2 == iid), ] <- PRS_info_matrix[which(fam$V2 == iid), ] - annotated_info$mean 
                PRS_info <- as.data.frame(t(PRS_info_matrix))
                colnames(PRS_info) <- fam$V2
                rownames(PRS_info) <- NULL
                PRS_info$Gene <- rep(annotated_info$symbol, length.out = nrow(PRS_info))
                end_time <- Sys.time()
                print(paste0("Time it takes for do some calculation for individual is ", end_time - start_time, " seconds"))


                list_of_objects <- list(annotated_info = annotated_info, full_prs = full_PRS, PRS_info = PRS_info, snp_info = snp, important_genes = important_genes)
                saveRDS(list_of_objects, "./output/indi_data.RDS")
           
                start_time <- Sys.time()
                which_symbol_is_important_individual(iid, indi_file, snp, full_PRS, PRS_info, annotated_info, important_genes,'./static/data/snp_level_manhattan_plot.json')
                end_time <- Sys.time()
                print(paste0("Time it takes for post-processing for individual is ", end_time - start_time, " seconds"))
                full_PRS <- full_PRS[-nrow(full_PRS), ]
                end_time_T <- Sys.time()
                time_taken_T <- end_time_T - start_time_T
                print(paste0("Total time it takes is ", time_taken_T, " seconds"))
            
                # Stop redirecting output
                sink()
        
            ''')

            return redirect(url_for("result_individual_2"))
        
        return render_template('result_2.html')
    return render_template('result_2.html')




@app.route('/run_1/result_individual' , methods=('GET', 'POST'))
def result_individual():
    return render_template('result_individual.html')

@app.route('/run_2/result_individual_2')
def result_individual_2():
    return render_template('result_individual_2.html')


@app.route('/delete_files', methods=['POST'])
def delete_files():
    delete_all_files()
    return jsonify({"message": "All files deleted from static/data and static/data/debug folders."})

# Register the delete_all_files function to be called on exit
atexit.register(delete_all_files)

if __name__=='__main__':
    app.run(debug=True, port=5000)
