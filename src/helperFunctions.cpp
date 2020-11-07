#include "helperFunctions.h"

void readMetaInfo(std::string file_path, std::vector<std::string> &individuals, std::vector<std::string> &site){
    std::cout << "Reading MetaInfo from the VCF file ...\n";
    std::string line = "";
    std::ifstream inFile(file_path);

    while (getline(inFile, line))
    {
        // READING THE INDIVIDUAL IDs FROM THE HEADER
        // more clear implementation with boost <TODO>
        if (line.find("#") != std::string::npos){
            if(line.find("#CHROM") != std::string::npos){
                std::istringstream is1(line);
                std::string wd;
                int32_t pos = 0;
                while( is1 >> wd){
                    if(pos < 9){

                        pos += 1;
                        continue;
                    }
                    else{
                        // to differentiate between the two haplotypes of the same individual
                        std::string temp0 = wd + "0";
                        individuals.push_back(temp0);
                        temp0 = wd + "1";
                        individuals.push_back(temp0);
                        pos += 1;
                    }
                }
            }
            continue;

        } else { 
            std::istringstream is(line);
            std::string wrd;
            int32_t pos = 0; //move through columns of each line of vcf file
            while(is >> wrd){
                if (pos == 1){
                    site.push_back(wrd);
                }
                ++pos;
            }
        }
    }// end reading vcf file

    std::cout << "Completed Reading the Meta-data!" << "\n";
    std::cout << "Total # of Individuals = " << individuals.size() << "\n";
    std::cout << "Total # of sites = " << site.size() << "\n";

    inFile.close();
}


std::vector<char> readGenotype(std::ifstream &inFile){
    std::string line = "";
    std::vector<char> X;
    while(getline(inFile, line)){
        if (line.find("#") != std::string::npos){
            continue;
        } else {
            std::istringstream is(line);
            std::string wrd;
            int32_t pos = 0;

            while(is >> wrd){
                if (pos >= 9){
                    X.push_back(wrd[0]);
                    X.push_back(wrd[2]);
                }
                ++pos;
            }
            break;
        }
    }
    return X;
}
