#ifndef DOUBLEPBWT_H
#define DOUBLEPBWT_H

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>

void doublePBWT(std::string vcfFile, std::vector<std::string> &individuals, std::vector<std::string> &site,
                double L, std::string genMapPath, int32_t overlap);

void startScan(std::string vcfFile, std::vector<std::string> &individuals, std::vector<std::string> &site,
               double L, std::string genMapPath);

void startScanOverlap(std::string vcfFile, std::vector<std::string> &individuals, std::vector<std::string> &site,
                      double L, std::string genMapPath, int32_t overlap);

// Strict Boundary
void leadingPBWT(std::vector<char> &X_lead, int32_t k_lag, std::vector<int32_t> &ppa_lead,
                 std::vector<int32_t> &div_lead, int32_t k_lead, int32_t N, std::vector<int32_t> &block_info, 
                 std::vector<int32_t> &group_info);

// Overlap 
void leadingPBWT(std::vector<char> &X_lead, int32_t k_lag, std::vector<int32_t> &ppa_lead,
                 std::vector<int32_t> &div_lead, int32_t k_lead, int32_t N, std::vector<int32_t> &block_info, 
                 std::vector<int32_t> &group_info, int32_t overlap);

void laggingAlgo3(std::vector<char> &X_lag, double L, std::vector<int32_t> &ppa_lag,
                  std::vector<int32_t> &div_lag, int32_t k_lag, int32_t N, std::vector<int32_t> &block_info,
                  std::vector<int32_t> &group_info, std::vector<std::string> &individuals,
                  std::vector<std::string> &site, std::map<int32_t, double> &genDist, int32_t k_lead, int32_t overlap, 
                  std::vector<int32_t> &div_lead, std::vector<int32_t> &ppa_lead);

// Report Alternating Matches for Strict Boundary
void findAltMatchesLessInfo(std::vector<int32_t> &ma, std::vector<int32_t> &mb,
                            std::vector<int32_t> &block_info, std::vector<int32_t> &group_info, int32_t k_lag,
                            std::vector<std::string> &individuals, std::vector<std::string> &site,
                            int32_t k_lead);
// Overlap
void findAltMatchesLessInfo(std::vector<int32_t> &ma, std::vector<int32_t> &mb,
                            std::vector<int32_t> &block_info, std::vector<int32_t> &group_info,
                            int32_t k_lag, std::vector<std::string> &individuals, std::vector<std::string> &site,
                            int32_t k_lead, int32_t overlap, std::vector<int32_t> &div_lead,
                            std::vector<int32_t> &ppa_lead);
// For strict Boundary
void updateBlockGroup(std::vector<int32_t> &div_lead, std::vector<int32_t> &ppa_lead,
                      int32_t N, std::vector<int32_t> &local_block_info,
                      std::vector<int32_t> &local_group_info, std::vector<std::string> &individuals, 
                      std::vector<std::string> &site, int32_t k_lag, int32_t k_lead);

// For overlap case
void updateBlockGroup(std::vector<int32_t> &div_lead, std::vector<int32_t> &ppa_lead,
                      int32_t N, std::vector<int32_t> &local_block_info, std::vector<int32_t> &local_group_info,
                      std::vector<std::string> &individuals, std::vector<std::string> &site,
                      int32_t k_lag, int32_t k_lead, int32_t overlap);

std::map<int32_t, double> getGenDist(std::string genMapPath);

int32_t findStartPosMatch(int32_t indv1, int32_t indv2, const std::vector<int32_t> &div,
                          const std::vector<int32_t> &ppa);
#endif
