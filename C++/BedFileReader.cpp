// header file
#include <fstream>  
#include <iostream> 
#include <sstream>
#include <stdio.h>
#include <bitset>
#include <vector>
#include <cstring>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <Rcpp.h>
#include "BedFileReader.h"

using namespace std;

/// @brief 
/// @param famName 
/// @param bimName 
/// @param bedName  




BedFileReader::BedFileReader(string famName, string bimName, string bedName){

    /* ---------- initialization ----------*/
    this->famName_temp = famName;
    this->bimName_temp = bimName;
    this->bedName_temp = bedName;
    this->m_line_counter = -1; // total number of individual 
    this->m_num_of_snps = -1; // total number of snps


    /* ----------iid_indexing ----------*/

    this->fam.open(this->famName_temp);

    // error messege
    if (!this->fam.is_open()){
    error.open("reading_file.log", ios_base::app);
    error << "FAM_FILE_READ_ERROR: cannot open the file\n";
    error.close();
    } 

    size_t curLine_2 = 0;
    string line_2;

    while(getline(this->fam, line_2)) { 
        curLine_2++;
        // cout << curLine_2 << endl;
        std::istringstream iss(line_2);

        int column3;
        string column4;

        iss >> column3 >> column4; 
        this->iid_index.insert({column4, curLine_2});
    }
    
    cout << "Done with iid indexing" << endl;

    this->fam.close(); 

    this->m_line_counter = curLine_2;
    cout << "Number of individuals: " << this->m_line_counter << endl;
    this->m_size_of_esi = (this->m_line_counter + 3)/4;
    // cout << "Number of bytes: " << this->m_size_of_esi << endl;


    /* ----------snp_indexing ----------*/
    this->bim.open(bimName);

    // error messege
    if (!this->bim.is_open()){
    error.open("reading_file.log", ios_base::app);
    error << "BIM_FILE_READ_ERROR: cannot open the file\n";
    error.close();
    } 

    size_t curLine = 0;
    string line;

    while(getline(this->bim, line)) { 
        curLine++;
        std::istringstream iss(line);

        int column1;
        string column2;

        iss >> column1 >> column2; 
        this->snp_index.insert({column2, curLine});
    }

    cout << "Done with snp indexing" << endl;
    this->bim.close(); 
    this->m_num_of_snps = curLine; 

    cout << "Number of snps: " << this->m_num_of_snps << endl;

    /* ---------- bed file reading ----------*/
        
    this->bed.open(bedName, ios::binary); // read it as binary

    //error messege
    if (!this->bed){
        error.open("reading_file.log", ios_base::app); 
        error << "CANT_OPEN_BED_FILE4READ\n";
        error.close();

    }

    char tmpbuf[3]; 
    memset(tmpbuf, '\0', sizeof(tmpbuf));
    this->bed.read(tmpbuf, sizeof(tmpbuf)); // three first bytes - permanent in bed file

    

}


BedFileReader::~BedFileReader() {
    close_bed();
}


void BedFileReader::close_bed() {
    if (fam.is_open()) fam.close();
    if (bim.is_open()) bim.close();
    if (bed.is_open()) bed.close();
}

// 리딩할때 고친거
vector<int> BedFileReader::readOneSnp_Whole(int snpIndex){


    this->snpIndex = snpIndex;
    /* ---------- initialization ----------*/
    int bits_val[MY_CHAR_BIT];
    vector<int> temp_snp_info0_4_snp(m_line_counter, 0);
    vector<int> temp_snp_info1_4_snp(m_line_counter, 0); 
    vector<char> encoded_snp_info_4_snp(m_size_of_esi, 0); 
    vector<char> buff_4_snp(m_size_of_esi, 0); 
    
    size_t individuals_counter = 0;
    size_t pos_4_snp = (snpIndex-1) * this->m_size_of_esi + 3;
    size_t pos_cur_4_snp = this->bed.tellg();

    // error messege
    if (!this->bed.is_open()){
    error.open("reading_file.log", ios_base::app);
    error << "1. BED_FILE_READ_ERROR: cannot open the file\n";
    error.close();
    } 

    if(!pos_cur_4_snp){
        error.open("reading_file.log", ios_base::app); 
        error << "2. BED_FILE_READ_ERROR: cannot get the position\n";
        error.close();
    }

    if(!this->bed.seekg(pos_4_snp, ios::beg)){
        error.open("reading_file.log", ios_base::app); 
        error << "3. BED_FILE_READ_ERROR: cannot get seek the position\n";
        error.close();
    }

    this->bed.read((char*) &buff_4_snp[0],this->m_size_of_esi);

    for (size_t i = 0; i < this->m_size_of_esi; i++) {	
    //===============================================================					
    //=== This part converts Byte "buff_4_snp[i]" to bits values "bits_val"
    //=== for example byte buff_4_snp[0] = "w" ->  bits_val = 11101011

        memset(bits_val, 0, sizeof(bits_val));

        int k = MY_CHAR_BIT;  //8
            while (k > 0){
                -- k;
                bits_val[k] = (buff_4_snp[i]&(1 << k) ? 1 : 0);
            }
    
        //=== here interpret Bit information "bits_val" to snps and count it - decode it
        decode_byte_4_snp(bits_val, &individuals_counter, &temp_snp_info0_4_snp[0], &temp_snp_info1_4_snp[0]);
        

    }

  return temp_snp_info0_4_snp;
}

vector<int> BedFileReader::readOneSnp(std::ifstream &bedFile, int snpIndex) {
    /* ---------- initialization ----------*/
    int bits_val[MY_CHAR_BIT];
    vector<int> temp_snp_info0_4_snp(m_line_counter, 0); // plink alternative allele
    vector<int> temp_snp_info1_4_snp(m_line_counter, 0); 
    vector<char> buff_4_snp(m_size_of_esi, 0);  // Buffer for SNP data

    // Calculate the position of the SNP in the .bed file
    size_t pos_4_snp = (snpIndex - 1) * this->m_size_of_esi + 3;  // +3 skips the first three magic bytes

    // Seek to the SNP position in the .bed file
    bedFile.seekg(pos_4_snp, ios::beg);
    bedFile.read(reinterpret_cast<char*>(&buff_4_snp[0]), this->m_size_of_esi);

    //고친것
    size_t individuals_counter = 0;

    // Decode SNP genotype data from buffer
    for (size_t i = 0; i < this->m_size_of_esi; i++) {
        memset(bits_val, 0, sizeof(bits_val));

        // Convert byte to bits (similar to before)
        for (int k = 0; k < MY_CHAR_BIT; ++k) {
            bits_val[k] = (buff_4_snp[i] & (1 << k)) ? 1 : 0;
        }
        // 고친것
        // size_t individuals_counter = i;  
        // Decode the bits to SNP genotypes
        decode_byte_4_snp(bits_val, &individuals_counter, &temp_snp_info0_4_snp[0], &temp_snp_info1_4_snp[0]);  // Fix: 4 arguments now
    
    }

    return temp_snp_info0_4_snp;
}



vector<int> BedFileReader::readOneIndi(string indiName) {
    auto indi_itr = this->iid_index.find(indiName);
    if (indi_itr == this->iid_index.end()) {
        cerr << "Individual ID not found" << endl;
        return vector<int>();
    }

    this->indiIndex = indi_itr->second;

    /* ---------- initialization ----------*/
    int bits_val[MY_CHAR_BIT];  // 8 bits per byte
    vector<int> temp_snp_info0_4_indi(m_num_of_snps, 0);  // SNP data for this individual
    vector<int> temp_snp_info1_4_indi(m_num_of_snps, 0); 
    vector<char> buff_4_indi(1, 0);  // Read 1 byte (since we're reading individuals in chunks of 4 per byte)

    size_t snps_counter = 0;

    // Error message if BED file isn't open
    if (!this->bed.is_open()) {
        error.open("reading_file.log", ios_base::app);
        error << "1. BED_FILE_READ_ERROR: cannot open the file\n";
        error.close();
    }

    for (size_t j = 0; j < m_num_of_snps; j++) {
        // Calculate the position of the byte that contains the genotype for the individual
        size_t pos_4_indi = floor((this->indiIndex - 1) / 4) + 3 + m_size_of_esi * j;

        // Seek to the correct byte in the .bed file
        this->bed.seekg(pos_4_indi, ios::beg);  // Ensure the file pointer moves correctly
        this->bed.read((char*)&buff_4_indi[0], 1);  // Read one byte

        // Clear the bits array
        memset(bits_val, 0, sizeof(bits_val));

        // Decode the byte into bits
        int k = MY_CHAR_BIT;  // 8 bits per byte
        while (k > 0) {
            --k;
            bits_val[k] = (buff_4_indi[0] & (1 << k)) ? 1 : 0;  // Extract bits
        }

        // Decode the SNP for this individual
        decode_byte_4_indi(bits_val, indiIndex, &snps_counter, &temp_snp_info0_4_indi[0], &temp_snp_info1_4_indi[0]);
    }

    // Reset the file pointer to the start of the .bed file for future operations
    this->bed.clear();
    this->bed.seekg(0, ios::beg);

    return temp_snp_info0_4_indi;
}

   

vector<float> BedFileReader::calculateWholePRS(
    vector<string>& snpList, 
    vector<float>& betaList, 
    const vector<int>& flipList) 
{
    // Check that all input vectors are of the same size
    if (snpList.size() != betaList.size() || snpList.size() != flipList.size()) {
        throw invalid_argument("snpList, betaList, and flipList must be of the same length.");
    }

    // Initialize PRS vector with zeros for all individuals
    vector<float> PRS(m_line_counter, 0.0f);  

    for (size_t i = 0; i < snpList.size(); i++) {
        // Find SNP index using snp_index map
        auto it = this->snp_index.find(snpList[i]);
        if (it == this->snp_index.end()) {
            // SNP not found in index map; skip or handle as needed
            cerr << "Warning: SNP " << snpList[i] << " not found in index map. Skipping.\n";
            continue;
        }
        int snpIndex = it->second;

        // Read the SNP genotype data for all individuals
        vector<int> oneSnp = readOneSnp_Whole(snpIndex);

        // Check that oneSnp has the correct number of individuals
        if (oneSnp.size() != m_line_counter) {
            cerr << "Error: Genotype data size mismatch for SNP " << snpList[i] << ".\n";
            continue;
        }

        // Determine if flipping is required for this SNP
        bool shouldFlip = (flipList[i] == 1);

        // Update PRS for each individual
        for (size_t j = 0; j < oneSnp.size(); j++) {
            int genotype = oneSnp[j];

            // Handle missing data (assuming 9 represents missing)
            if (genotype == 9) {
                continue;  // Skip missing genotypes
            }

            // Apply flipping if required
            if (shouldFlip) {
                genotype = 2 - genotype;
            }

            // Accumulate PRS
            PRS[j] += static_cast<float>(genotype) * betaList[i];
        }

        // Reset the file pointer between SNP reads
        this->bed.clear();
        this->bed.seekg(0, ios::beg);
    }

    return PRS;
}






std::vector<float> BedFileReader::calculatePRS(
    const std::vector<std::string>& snpList, 
    const std::vector<float>& betaList, 
    const std::vector<int>& flipList,
    const std::vector<float>& AFList) 
{
    // Check that all input vectors are of the same size
    if (snpList.size() != betaList.size() || 
        snpList.size() != flipList.size() ||
        snpList.size() != AFList.size()) {  // Check AFList size
        throw std::invalid_argument("snpList, betaList, flipList, and AFList must all be of the same length.");
    }

    // Initialize PRS vector with zeros for all individuals
    std::vector<float> PRS(m_line_counter, 0.0f);  

    // Open the .bed file separately for each process/thread
    std::ifstream bedFile(this->bedName_temp, ios::binary);
    if (!bedFile.is_open()) {
        std::ofstream error("reading_file.log", ios_base::app);
        error << "BED_FILE_READ_ERROR: cannot open the .bed file\n";
        error.close();
        return PRS;  // Return empty PRS vector if the file can't be opened
    }

    // For each SNP in snpList, calculate PRS
    for (size_t i = 0; i < snpList.size(); i++) {
        // Find SNP index
        auto snpIndexItr = this->snp_index.find(snpList[i]);
        if (snpIndexItr == this->snp_index.end()) {
            std::cerr << "Warning: SNP " << snpList[i] << " not found in index map. Skipping.\n";
            continue;  // Skip if SNP is not found
        }

        int snpIndex = snpIndexItr->second;

        // Read the SNP genotype data for all individuals
        std::vector<int> oneSnp = readOneSnp(bedFile, snpIndex);  // Pass the bed file descriptor to readOneSnp

        // Check that oneSnp has the correct number of individuals
        if (oneSnp.size() != m_line_counter) {
            std::cerr << "Error: Genotype data size mismatch for SNP " << snpList[i] << ".\n";
            continue;
        }

         // Retrieve allele frequency for this SNP
        float alleleFreq = AFList[i];

        // Determine if flipping is required for this SNP
        bool shouldFlip = (flipList[i] == 1);

        // Update PRS for each individual
        for (size_t j = 0; j < oneSnp.size(); j++) {
            int genotype = oneSnp[j];

            // Handle missing data (assuming 9 represents missing)
            if (genotype == 9) {
                continue;  // Skip missing genotypes
            }

            float adjustedGenotype;

            if (shouldFlip) {
                // Apply flipping: 2 - genotype
                adjustedGenotype = static_cast<float>(2 - genotype);
            }
            else {
                // No flipping
                adjustedGenotype = static_cast<float>(genotype);
            }

            // Subtract allele frequency
            adjustedGenotype -= alleleFreq;

            // Accumulate PRS
            PRS[j] += adjustedGenotype * betaList[i];

            //Debugging log
            // std::cerr << "SNP " << snpList[i] << ", individual " << j 
            //           << ": genotype=" << genotype 
            //           << ", adjustedGenotype=" << adjustedGenotype 
            //           << ", PRS=" << PRS[j] << "\n";

            //  std::cerr << "SNP " << snpList[i] << ": genotype=" << genotype << endl;
                     
        
            }
        }

    bedFile.close();  // Close the file when done
    return PRS;
}


void BedFileReader::decode_byte_4_snp(int* bits_val, size_t * individuals_counter, int* temp_snp_info0, int* temp_snp_info1){
    for (int i = 0; i < 4; ++i)
    {
        if(*individuals_counter >= m_line_counter)
            return; // Prevent overflow

        if(bits_val[i*2] == 0 && bits_val[i*2+1] == 0)
        {
            temp_snp_info0[*individuals_counter] = 2; // Homozygote 1
            temp_snp_info1[*individuals_counter] = 2;
        }
        else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 1)
        {
            temp_snp_info0[*individuals_counter] = 0; // Homozygote 2
            temp_snp_info1[*individuals_counter] = 0;
        }
        else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 0)
        {
            temp_snp_info0[*individuals_counter] = 9; // Missing genotype
            temp_snp_info1[*individuals_counter] = 9;
        }
        else if(bits_val[i*2] == 0 && bits_val[i*2+1] == 1)
        {
            temp_snp_info0[*individuals_counter] = 1; // Heterozygote
            temp_snp_info1[*individuals_counter] = 1;
        }

        *individuals_counter += 1; // Increment after processing each genotype
    }
}


void BedFileReader::decode_byte_4_indi(int* bits_val, int indiIndex, size_t* snps_counter, int* temp_snp_info0, int* temp_snp_info1) {
    // Calculate which 2-bit pair (out of the 4 in the byte) corresponds to this individual
    int i = (indiIndex - 1) % 4;  // 0-based index

    // Decode the genotype for this individual
    if (bits_val[i * 2] == 0 && bits_val[i * 2 + 1] == 0) {
        *snps_counter += 1;
        if (*snps_counter > m_num_of_snps)
            return;
        temp_snp_info0[*snps_counter - 1] = 2;  // Homozygous for A2 allele
        temp_snp_info1[*snps_counter - 1] = 0;
    }
    else if (bits_val[i * 2] == 1 && bits_val[i * 2 + 1] == 1) {
        *snps_counter += 1;
        if (*snps_counter > m_num_of_snps)
            return;
        temp_snp_info0[*snps_counter - 1] = 0;  // Homozygous for A1 allele
        temp_snp_info1[*snps_counter - 1] = 2;
    }
    else if (bits_val[i * 2] == 1 && bits_val[i * 2 + 1] == 0) {
        *snps_counter += 1;
        if (*snps_counter > m_num_of_snps)
            return;
        temp_snp_info0[*snps_counter - 1] = 9;  // Missing value
        temp_snp_info1[*snps_counter - 1] = 9;
    }
    else if (bits_val[i * 2] == 0 && bits_val[i * 2 + 1] == 1) {
        *snps_counter += 1;
        if (*snps_counter > m_num_of_snps)
            return;
        temp_snp_info0[*snps_counter - 1] = 1;  // Heterozygous A1/A2
        temp_snp_info1[*snps_counter - 1] = 1;
    }
}



   
// [[Rcpp::export]]
/// create an external pointer to a Uniform object
RcppExport SEXP BedFileReader__new(SEXP famName_, SEXP bimName_, SEXP bedName_) {
  // convert inputs to appropriate C++ types
  string famName = Rcpp::as<string>(famName_);
  string bimName = Rcpp::as<string>(bimName_);
  string bedName = Rcpp::as<string>(bedName_);

  // create a pointer to an Uniform object and wrap it as an external pointer
  Rcpp::XPtr<BedFileReader> reader( new BedFileReader(famName, bimName, bedName), true );

  // return the external pointer to the R side
  return reader;
}


RcppExport SEXP BedFileReader__readOneSnp_Whole(SEXP xp, SEXP snpIndex_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    int snpIndex = Rcpp::as<int>(snpIndex_);

    vector<int> oneSnp = reader->readOneSnp_Whole(snpIndex);
    
    return Rcpp::wrap(oneSnp);
}


RcppExport SEXP BedFileReader__readOneIndi(SEXP xp, SEXP indiName_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    string indiName = Rcpp::as<string>(indiName_);

    vector<int> oneIndi = reader->readOneIndi(indiName);
    
    return Rcpp::wrap(oneIndi);
}

RcppExport SEXP BedFileReader__calculatePRS(SEXP xp, SEXP snpList_, SEXP betaList_, SEXP flipList_, SEXP AFList_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    vector<string> snpList = Rcpp::as< vector<string> >(snpList_);
    vector<float> betaList = Rcpp::as< std::vector<float> >(betaList_);
    vector<int> flipList = Rcpp::as< std::vector<int> >(flipList_);
    vector<float> AFList = Rcpp::as< std::vector<float> >(AFList_);  // Extract AFList

     // Ensure that the sizes of snpList, betaList, flipList, and AFList match
    if (snpList.size() != betaList.size() || 
        snpList.size() != flipList.size() ||
        snpList.size() != AFList.size()) {
        Rcpp::stop("snpList, betaList, flipList, and AFList must all be of the same length.");
    }
    // Calculate PRS using the updated function
     std::vector<float> PRS = reader->calculatePRS(snpList, betaList, flipList, AFList);


    // Return PRS as NumericVector to R
    return Rcpp::wrap(PRS);
}

RcppExport SEXP BedFileReader__calculateWholePRS(SEXP xp, SEXP snpList_, SEXP betaList_, SEXP flipList_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    vector<string> snpList = Rcpp::as< vector<string> >(snpList_);
    vector<float> betaList = Rcpp::as< std::vector<float> >(betaList_);
    vector<int> flipList = Rcpp::as< std::vector<int> >(flipList_);

    // Ensure that the sizes of snpList, betaList, and flipList match
    if (snpList.size() != betaList.size() || snpList.size() != flipList.size()) {
        Rcpp::stop("snpList, betaList, and flipList must be of the same length.");
    }

    vector<float> PRS = reader->calculateWholePRS(snpList, betaList, flipList);

    return Rcpp::wrap(PRS);
}


