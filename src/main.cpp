#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "doublePBWT.h"
#include "helperFunctions.h"

// ************** Main ************ //
int main(int argc, char *argv[])
{
    std::ios_base::sync_with_stdio(false);
    std::string vcfFile = "";
    std::string genMapPath = "";
    double L = 0;
    int32_t overlap = 0;
    std::vector<std::string> individuals, site;

    if (argc == 5)
    {
        vcfFile = std::string(argv[1]);
        genMapPath = std::string(argv[2]);
        L = std::stod(argv[3]);
        overlap = std::stoi(argv[4]);
    }
    else
    {
        std::cout << "Usage: \n";
        std::cout << "./d2pbwt_overlap_genetic <vcf-file-path> <gen-map-filepath> <gendist-L-threshold> <max-overlap-in-sites>\n";
        exit(1);
    }

    // populate individuals and site vector
    readMetaInfo(vcfFile, individuals, site);

    std::cout << "gendist L = " << L << "cM\n";
    std::cout << "overlap = " << overlap << "sites\n";

    doublePBWT(vcfFile, individuals, site, L, genMapPath, overlap);
    return 0;
}
