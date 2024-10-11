// header file
#include <fstream>  
#include <iostream> 
#include <sstream>
#include <stdio.h>
#include <bitset>
#include <vector>
#include <cstring>
#include <chrono>
#include <cmath>
#include <Rcpp.h>
#include "BedFileReader_indi.h"

#ifndef _GLIBCXX_USE_C99_TIMESPEC
#define _GLIBCXX_USE_C99_TIMESPEC 1
#endif

using namespace std;

/// @brief 
/// @param famName 
/// @param bimName 
/// @param bedName  


BedFileReader_indi::BedFileReader_indi(string famName, string bimName, string bedName) {
    famName_temp = famName;
    bimName_temp = bimName;
    bedName_temp = bedName;

    fam.open(famName);
    bim.open(bimName);
    bed.open(bedName, ios::in | ios::binary);

    // Check if files are opened successfully
    if (!fam.is_open()) {
        cerr << "FAM_FILE_READ_ERROR: cannot open the file" << endl;
    }
    if (!bim.is_open()) {
        cerr << "BIM_FILE_READ_ERROR: cannot open the file" << endl;
    }
    if (!bed.is_open()) {
        cerr << "BED_FILE_READ_ERROR: cannot open the file" << endl;
    }
}

BedFileReader_indi::~BedFileReader_indi() {
    close_bed();
}

void BedFileReader_indi::close_bed() {
    if (fam.is_open()) fam.close();
    if (bim.is_open()) bim.close();
    if (bed.is_open()) bed.close();
}

void BedFileReader_indi::initialize_indices() {
    string line;
    int index = 0;

    // Reopen FAM file to build iid_index
    fam.clear();
    fam.seekg(0, ios::beg);
    while (getline(fam, line)) {
        istringstream iss(line);
        string fid, iid;
        iss >> fid >> iid;  // Read family ID and individual ID
        iid_index[iid] = index++;
    }
    m_line_counter = index;

    // Debugging output: Print parsed individual IDs and indices
    // cout << "Parsed individual IDs and indices:" << endl;
    // for (const auto& pair : iid_index) {
    //     cout << pair.first << " -> " << pair.second << endl;
    // }

    // Reopen BIM file to build snp_index
    index = 0;
    bim.clear();
    bim.seekg(0, ios::beg);
    while (getline(bim, line)) {
        istringstream iss(line);
        string snp;
        for (int i = 0; i < 2; ++i) iss >> snp;  // Skip first two columns
        iss >> snp;  // Read SNP ID
        snp_index[snp] = index++;
    }
    m_num_of_snps = index;

    m_size_of_esi = (m_line_counter + 3) / 4;

    // Debugging output: Verify the number of SNPs and size of ESI
    cout << "Number of SNPs: " << m_num_of_snps << endl;
    cout << "Number of Individuals: " << m_line_counter << endl;
}



vector<int> BedFileReader_indi::readOneSnp(int snpIndex){
     if (iid_index.empty() || snp_index.empty()) {
        initialize_indices();
    }


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

vector<int> BedFileReader_indi::readOneIndi(string iid) {
    if (iid_index.empty() || snp_index.empty()) {
        initialize_indices();
    }

    auto indi_itr = this->iid_index.find(iid);
    if (indi_itr == this->iid_index.end()) {
        cerr << "Individual ID not found: " << iid << endl;
        return vector<int>();
    }

    this->indiIndex = indi_itr->second;

    /* ---------- initialization ----------*/
    vector<int> snp_data(m_num_of_snps, -1);
    int bits_val[8];  // Define bits_val with size 8
    vector<char> buff_4_indi(m_size_of_esi, 0); 

    size_t pos_4_indi_base = floor(this->indiIndex / 4) + 3;

    // Error message handling
    if (!this->bed.is_open()) {
        cerr << "1. BED_FILE_READ_ERROR: cannot open the file" << endl;
        return vector<int>();
    }

    for (size_t snp = 0; snp < m_num_of_snps; ++snp) {
        size_t pos_4_indi = pos_4_indi_base + m_size_of_esi * snp;
        this->bed.seekg(pos_4_indi, ios::beg);
        this->bed.read((char*)&buff_4_indi[0], 1);

        if (!this->bed) {
            cerr << "4. BED_FILE_READ_ERROR: reading error" << endl;
            return vector<int>();
        }

        unsigned char byte = buff_4_indi[0];
        int bit_pos = (this->indiIndex % 4) * 2;  // Each individual is represented by 2 bits in the byte
        int genotype = (byte >> bit_pos) & 0b11;  // Extract the two bits for the genotype

        switch (genotype) {
            case 0b00:  // Homozygous major
                snp_data[snp] = 2;
                break;
            case 0b11:  // Homozygous minor
                snp_data[snp] = 0;
                break;
            case 0b10:  // Heterozygous
                snp_data[snp] = 1;
                break;
            case 0b01:  // Missing genotype
            default:
                snp_data[snp] = -1;
                break;
        }
    }




    return snp_data;
}



// void BedFileReader::readAllSnp(string fileName){

//     chrono::steady_clock::time_point begin = chrono::steady_clock::now();
//     ofstream matrix;
//     matrix.open(fileName, ios_base::app);
//     for (int i= 1; i <= this->m_num_of_snps; i++){
//         vector<int> a0 = readOneSnp(i); 
//         for (int j= 0; j < this->m_line_counter; j++){
//             matrix << a0[j] << " ";
//         }
//         matrix << "\n";
        
//     }
//     matrix.close();
//     chrono::steady_clock::time_point end = chrono::steady_clock::now();
//     cout << "Time take to read all snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;

// }

vector<float> BedFileReader_indi::calculatePRS(vector<string> snpList, vector<float> betaList){

    vector<int> oneSnp; 
    vector<float> PRS(m_line_counter, 0);


    for (int i = 0; i < snpList.size(); i++){
        //snpIndex = findSnpIndex(snpList[i]);
        
        snpIndex = this->snp_index.find(snpList[i])->second;
        oneSnp = readOneSnp(snpIndex);
        for (int j = 0; j < oneSnp.size(); j++){
            PRS[j] += oneSnp[j] * betaList[i];
        }
    }
    return PRS;
}


void BedFileReader_indi::decode_byte_4_snp(int* bits_val, size_t * individuals_counter, int* temp_snp_info0, int* temp_snp_info1){
	//int flag = 0;
	for (int i = 0; i < 4; ++i)
	{
		if(bits_val[i*2] == 0 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 2;
			temp_snp_info1[*individuals_counter - 1] = 0;
			//flag = 1 ;//Homozegote 1 for example GG      write 20 ; 00 - will point to [0] letter   + 2 to [0]

		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 0;
			temp_snp_info1[*individuals_counter - 1] = 2;
			//flag = 2 ;//Homozegote 2 for example AA      write 02 ; 11 - will point to [1] letter   + 2 to [1]
		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;		
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 9;
			temp_snp_info1[*individuals_counter - 1] = 9;
			//flag = 3 ; //Missing value                   nothing to add - write 99 ;
		}
		else if(bits_val[i*2] == 0 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 1;
			temp_snp_info1[*individuals_counter - 1] = 1;

			//flag = 4 ; //Heterozegote for example AG  or GA     write 11 ; 01 - will point to [0] and [1] letter   +1 +1 to [1]
		}
		// else
		// 	flag = 5 ; //Error
	}
}

 


   
// [[Rcpp::export]]
/// create an external pointer to a Uniform object
RcppExport SEXP BedFileReader_indi_new(SEXP famName_, SEXP bimName_, SEXP bedName_) {
  // convert inputs to appropriate C++ types
  string famName = Rcpp::as<string>(famName_);
  string bimName = Rcpp::as<string>(bimName_);
  string bedName = Rcpp::as<string>(bedName_);

  // create a pointer to an Uniform object and wrap it as an external pointer
  Rcpp::XPtr<BedFileReader_indi> reader( new BedFileReader_indi(famName, bimName, bedName), true );

  // return the external pointer to the R side
  return reader;
}


RcppExport SEXP BedFileReader_indi_readOneSnp(SEXP xp, SEXP snpIndex_){
    Rcpp::XPtr<BedFileReader_indi> reader(xp);

    int snpIndex = Rcpp::as<int>(snpIndex_);

    vector<int> oneSnp = reader->readOneSnp(snpIndex);
    
    return Rcpp::wrap(oneSnp);
}


RcppExport SEXP BedFileReader_indi_readOneIndi(SEXP xp, SEXP indiName_){
    Rcpp::XPtr<BedFileReader_indi> reader(xp);

    string indiName = Rcpp::as<string>(indiName_);

    vector<int> oneIndi = reader->readOneIndi(indiName);
    
    return Rcpp::wrap(oneIndi);
}

RcppExport SEXP BedFileReader_indi_calculatePRS(SEXP xp, SEXP snpList_, SEXP betaList_){
    Rcpp::XPtr<BedFileReader_indi> reader(xp);

    vector<string> snpList = Rcpp::as< vector<string> >(snpList_);
    vector<float> betaList = Rcpp::as< std::vector<float> >(betaList_);

    vector<float> PRS = reader->calculatePRS(snpList, betaList);

    return Rcpp::wrap(PRS);
}

