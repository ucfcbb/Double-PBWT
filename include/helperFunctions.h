#ifndef HELPER_H
#define HELPER_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

void readMetaInfo(std::string file_path, std::vector<std::string> &individuals, std::vector<std::string> &site);
std::vector<char> readGenotype(std::ifstream &inFile);

#endif
