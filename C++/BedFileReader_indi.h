#ifndef _BEDFILEREADER_H        
#define _BEDFILEREADER_H 

#include <fstream>  
#include <iostream> 
#include <map>
#include <vector>


using namespace std;


#define MY_CHAR_BIT 8


class BedFileReader_indi{

    public:
    /* ---------- Constructor ----------*/
    BedFileReader_indi(string famName, string bimName, string bedName);
    ~BedFileReader_indi();
    /* ---------- Member functions ----------*/
    vector<int> readOneSnp(int snpIndex);
    vector<int> readOneIndi(string indiName);
    // vector<int> readOneIndi(vector<string> snpList, string indiName);
    // void readAllSnp(string fileName);
    void close_bed();
    vector<float> calculatePRS(vector<string> snpList, vector<float> betaList);

    
    private:

    /* ---------- Data members ----------*/
    
    string famName_temp;
    string bimName_temp;
    string bedName_temp;
    
    // Streams for files
    ifstream fam;
    ifstream bim;
    ifstream bed;
    ofstream error;

    // Variables
    size_t m_line_counter; // Total number of individual 
    size_t m_size_of_esi; // Total number of bytes that hold all genotype info per one snp = one line length!!!
    size_t m_num_of_snps; // Total number of Snps
    map<string, int> snp_index;
    map<string, int> iid_index;
    int snpIndex;
    int indiIndex;

   

    /* ---------- Member functions ----------*/
    void initialize_indices();
    void decode_byte_4_snp(int* bits_val, size_t * individuals_counter, int* temp_snp_info0, int* temp_snp_info1);
    void decode_byte_4_indi(int* bits_val, int indiIndex ,size_t * snps_counter, int* temp_snp_info0, int* temp_snp_info1);

    };

    #endif //_BEDFILEREADER_H        