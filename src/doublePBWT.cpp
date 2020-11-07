#include "doublePBWT.h"
#include "helperFunctions.h"

#include <algorithm>

void doublePBWT(std::string vcfFile, std::vector<std::string> &individuals, std::vector<std::string> &site,
                double L, std::string genMapPath, int32_t overlap)
{
    if (overlap == 0)
    {
        // Juxtaposed Boundary case
        startScan(vcfFile, individuals, site, L, genMapPath);
    }
    else if (overlap > 0)
    {
        // Overlap Case
        startScanOverlap(vcfFile, individuals, site, L, genMapPath, overlap);
    }
    else
    {
        std::cout << "Overlap should be >= 0.\n";
        exit(1);
    }
}

void startScan(std::string vcfFile, std::vector<std::string> &individuals, std::vector<std::string> &site,
               double L, std::string genMapPath)
{
    std::ifstream inFile_lead(vcfFile);
    std::ifstream inFile_lag(vcfFile);

    int32_t k_lead = 0;
    int32_t k_lag = 0;
    const int32_t M = individuals.size();
    const int32_t N = site.size();
    std::vector<int32_t> ma_lag, mb_lag;
    std::vector<int32_t> ppa_lead, ppa_lag;
    std::vector<int32_t> div_lead(M, 0);
    std::vector<int32_t> div_lag(M, 0);

    std::vector<int32_t> block_info(M, 0);
    std::vector<int32_t> group_info(M, 0);
    std::vector<char> X_lead, X_lag;

    // intialize ppa
    for (int32_t i = 0; i < M; ++i)
    {
        ppa_lead.push_back(i);
        ppa_lag.push_back(i);
    }

    // populate GenDist() map
    std::map<int32_t, double> genDist = getGenDist(genMapPath);

    bool moveForwardOnce = true;
    while (k_lead <= N)
    {
        if (genDist[k_lag] < L)
        {
            // leading PBWT
            X_lead = readGenotype(inFile_lead);
            leadingPBWT(X_lead, k_lag, ppa_lead, div_lead, k_lead, N, block_info, group_info);
            k_lead += 1;

            for (auto i = 0; i < M; ++i)
            {
                group_info[i] = 0;
                block_info[i] = 0;
            }

            // run lagging PBWT
            X_lag = readGenotype(inFile_lag);
            laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, block_info, group_info, individuals,
                         site, genDist, k_lead, 0, div_lead, ppa_lead);
            // laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, block_info, group_info, individuals,
            //              site, genDist, k_lead);
            k_lag += 1;
            continue;
        }

        if (moveForwardOnce)
        {
            // leading PBWT
            X_lead = readGenotype(inFile_lead);
            leadingPBWT(X_lead, k_lag, ppa_lead, div_lead, k_lead, N, block_info,
                        group_info);
            k_lead += 1;
            for (auto i = 0; i < M; ++i)
            {
                group_info[i] = 0;
                block_info[i] = 0;
            }
            // run lagging PBWT
            X_lag = readGenotype(inFile_lag);
            laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, block_info, group_info, individuals,
                         site, genDist, k_lead, 0, div_lead, ppa_lead);
            k_lag += 1;
            moveForwardOnce = false;
            continue;
        }

        // leading PBWT
        X_lead = readGenotype(inFile_lead);
        leadingPBWT(X_lead, k_lag, ppa_lead, div_lead, k_lead, N, block_info, group_info);
        if (genDist[k_lead - 1] - genDist[k_lag] >= L)
        {
            X_lag = readGenotype(inFile_lag);
            laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, block_info, group_info, individuals,
                         site, genDist, k_lead, 0, div_lead, ppa_lead);
            // laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, block_info, group_info, individuals,
            //              site, genDist, k_lead);
            k_lag += 1;
        }
        else
        {
            // moving leading PBWT forward to satisfy length constraint
            k_lead += 1;
            for (auto i = 0; i < M; ++i)
            {
                group_info[i] = 0;
                block_info[i] = 0;
            }
            continue;
        }

        std::vector<int32_t> local_block_info(M, 0);
        std::vector<int32_t> local_group_info(M, 0);

        // Shifting lagging PBWT as long as the length constraint holds
        while (genDist[k_lead - 1] - genDist[k_lag] >= L)
        {
            updateBlockGroup(div_lead, ppa_lead, N, local_block_info, local_group_info,
                             individuals, site, k_lag, k_lead);
            // call lagging PBWT
            X_lag = readGenotype(inFile_lag);
            laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, local_block_info, local_group_info, individuals,
                         site, genDist, k_lead, 0, div_lead, ppa_lead);
            // laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, local_block_info,
            //              local_group_info, individuals, site, genDist, k_lead);
            k_lag += 1;
            for (auto i = 0; i < M; ++i)
            {
                local_group_info[i] = 0;
                local_block_info[i] = 0;
            }
        }

        // re-calibrate block_info and group_info for next k_lead
        // switch_info.clear();
        for (auto i = 0; i < M; ++i)
        {
            group_info[i] = 0;
            block_info[i] = 0;
        }
        k_lead += 1;
        X_lead.clear();
        X_lag.clear();
    }
}

void startScanOverlap(std::string vcfFile, std::vector<std::string> &individuals, std::vector<std::string> &site,
                      double L, std::string genMapPath, int32_t overlap)
{
    std::ifstream inFile_lead(vcfFile);
    std::ifstream inFile_lag(vcfFile);

    int32_t k_lead = 0;
    int32_t k_lag = 0;
    const int32_t M = individuals.size();
    const int32_t N = site.size();
    std::vector<int32_t> ma_lag, mb_lag;
    std::vector<int32_t> ppa_lead, ppa_lag;
    std::vector<int32_t> div_lead(M, 0);
    std::vector<int32_t> div_lag(M, 0);

    std::vector<int32_t> block_info(M, 0);
    std::vector<int32_t> group_info(M, 0);
    std::vector<char> X_lead, X_lag;

    // intialize ppa
    for (int32_t i = 0; i < M; ++i)
    {
        ppa_lead.push_back(i);
        ppa_lag.push_back(i);
    }

    // populate genDist map
    std::map<int32_t, double> genDist = getGenDist(genMapPath);
    bool moveForward = true;
    while (k_lead <= N)
    {
        if (genDist[k_lag] < L)
        {
            // leading PBWT
            X_lead = readGenotype(inFile_lead);
            leadingPBWT(X_lead, k_lag, ppa_lead, div_lead, k_lead, N, block_info, group_info, overlap);
            k_lead += 1;

            // recalibrate to avoid false alt-matches
            for (auto i = 0; i < M; ++i)
            {
                block_info[i] = 0;
            }

            // run lagging PBWT
            X_lag = readGenotype(inFile_lag);
            laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, block_info, group_info, individuals,
                         site, genDist, k_lead, overlap, div_lead, ppa_lead);
            k_lag += 1;
            continue;
        }

        if (moveForward)
        {
            // leading PBWT
            X_lead = readGenotype(inFile_lead);
            leadingPBWT(X_lead, k_lag, ppa_lead, div_lead, k_lead, N, block_info, group_info, overlap);
            k_lead += 1;
            for (auto i = 0; i < M; ++i)
            {
                block_info[i] = 0;
            }
            // run lagging PBWT
            X_lag = readGenotype(inFile_lag);
            laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, block_info, group_info, individuals,
                         site, genDist, k_lead, overlap, div_lead, ppa_lead);
            k_lag += 1;
            moveForward = false;
            continue;
        }

        // leading PBWT
        X_lead = readGenotype(inFile_lead);
        leadingPBWT(X_lead, k_lag, ppa_lead, div_lead, k_lead, N, block_info, group_info, overlap);
        if (genDist[k_lead - 1] - genDist[k_lag] >= L)
        {
            X_lag = readGenotype(inFile_lag);
            laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, block_info, group_info, individuals,
                         site, genDist, k_lead, overlap, div_lead, ppa_lead);
            k_lag += 1;
            X_lag.clear();
        }
        else
        {
            // moving only leading PBWT forward to satisfy gen. dist. constraint
            k_lead += 1;
            for (auto i = 0; i < M; ++i)
            {
                block_info[i] = 0;
            }
            continue;
        }

        // Shifting lagging PBWT as long as the length constraint holds
        std::vector<int32_t> local_block_info(M, 0);
        std::vector<int32_t> local_group_info(M, 0);
        while (genDist[k_lead - 1] - genDist[k_lag] >= L)
        {
            updateBlockGroup(div_lead, ppa_lead, N, local_block_info, local_group_info, individuals, site, k_lag, k_lead, overlap);

            // call lagging PBWT
            X_lag = readGenotype(inFile_lag);
            laggingAlgo3(X_lag, L, ppa_lag, div_lag, k_lag, N, local_block_info, local_group_info, individuals, site,
                         genDist, k_lead, overlap, div_lead, ppa_lead);
            k_lag += 1;
            for (auto i = 0; i < M; ++i)
            {
                local_block_info[i] = 0;
                local_group_info[i] = 0;
            }
        }

        // re-calibrate block_info and group_info for next k_lead
        for (auto i = 0; i < M; ++i)
        {
            block_info[i] = 0;
            group_info[i] = 0;
        }
        k_lead += 1;
        X_lead.clear();
        X_lag.clear();
    }
}

// Strict Boundary
void leadingPBWT(std::vector<char> &X_lead, int32_t k_lag, std::vector<int32_t> &ppa_lead,
                 std::vector<int32_t> &div_lead, int32_t k_lead, int32_t N, std::vector<int32_t> &block_info, std::vector<int32_t> &group_info)
{
    std::vector<int32_t> a, b, d, e;
    auto q = k_lead + 1;
    auto p = q;
    int32_t block_id = 1;
    bool find_block_start = true;
    std::vector<int32_t> curr_indiv;
    int32_t M = ppa_lead.size();

    // Boundary Case
    if (k_lead > N - 1)
    {
        // std::cout << "Boundary Case!! \n";
        for (int32_t index = div_lead.size() - 1; index > 0; --index)
        {
            if (find_block_start)
            {

                if (div_lead[index] < k_lead - (k_lead - k_lag))
                {
                    curr_indiv.push_back(ppa_lead[index]);
                }
                else if (div_lead[index] == k_lead - (k_lead - k_lag))
                {
                    curr_indiv.push_back(ppa_lead[index]);
                    for (auto indiv : curr_indiv)
                    {
                        block_info[indiv] = block_id;
                        group_info[indiv] = 1;
                    }

                    find_block_start = false;
                    curr_indiv.clear();
                }
                else
                {
                    curr_indiv.clear();
                    continue;
                }
            }
            else
            {
                // block boundary continues
                if (div_lead[index] < k_lead - (k_lead - k_lag))
                    curr_indiv.push_back(ppa_lead[index]);
                else
                {
                    // block boundary ends
                    curr_indiv.push_back(ppa_lead[index]);
                    for (auto indiv : curr_indiv)
                    {
                        block_info[indiv] = block_id;
                        group_info[indiv] = -1;
                    }
                    find_block_start = true;
                    curr_indiv.clear();
                    block_id += 1;
                }
            }
        }
        return;
    }

    for (auto i = 0; i < M; ++i)
    {
        auto div_index = M - 1 - i;
        if (div_index >= 0)
        {
            if (find_block_start)
            {
                // block boundary possible start
                if (div_lead[div_index] < k_lag)
                {
                    curr_indiv.push_back(ppa_lead[div_index]);
                }
                else if (div_lead[div_index] == k_lag)
                {
                    // VALID BLOCK
                    curr_indiv.push_back(ppa_lead[div_index]);
                    for (auto indiv : curr_indiv)
                    {
                        block_info[indiv] = block_id;
                        group_info[indiv] = 1;
                    }

                    find_block_start = false;
                    curr_indiv.clear();
                }
                else
                {
                    curr_indiv.clear();
                    // continue;
                }
            }
            else
            {
                // block boundary continues
                if (div_lead[div_index] < k_lag)
                    curr_indiv.push_back(ppa_lead[div_index]);
                else
                {
                    // block boundary ends
                    curr_indiv.push_back(ppa_lead[div_index]);
                    for (auto indiv : curr_indiv)
                    {
                        block_info[indiv] = block_id;
                        group_info[indiv] = -1;
                    }
                    find_block_start = true;
                    curr_indiv.clear();
                    block_id += 1;
                }
            }
        }

        // update ppa and divergence arrays
        auto index = ppa_lead[i];
        auto match_start = div_lead[i];
        auto allele = X_lead[index] - '0';

        if (match_start > p)
            p = match_start;
        if (match_start > q)
            q = match_start;

        if (allele == 0)
        {
            a.push_back(index);
            d.push_back(p);
            p = 0;
        }
        else
        {
            b.push_back(index);
            e.push_back(q);
            q = 0;
        }
    }

    if (k_lead < N)
    {
        ppa_lead = a;
        ppa_lead.insert(ppa_lead.end(), b.begin(), b.end());

        div_lead = d;
        div_lead.insert(div_lead.end(), e.begin(), e.end());
    }

    // std::cout << "End of Leading PBWT!\n";
}
// Overlap
void leadingPBWT(std::vector<char> &X_lead, int32_t k_lag, std::vector<int32_t> &ppa_lead,
                 std::vector<int32_t> &div_lead, int32_t k_lead, int32_t N,
                 std::vector<int32_t> &block_info, std::vector<int32_t> &group_info, int32_t overlap)
{
    std::vector<int32_t> a, b, d, e;
    auto q = k_lead + 1;
    auto p = q;
    int32_t block_id = 1, group_id = 1;
    std::vector<int32_t> temp_indiv;
    int32_t M = ppa_lead.size();
    bool start = true;
    // Boundary Case
    if (k_lead > N - 1)
    {
        // std::cout << "Boundary Case!! \n";
        for (int32_t index = div_lead.size() - 1; index >= 0; --index)
        {
            if (start)
            {
                // block boundary possible start
                if (div_lead[index] < k_lag - overlap)
                {
                    temp_indiv.push_back(ppa_lead[index]);
                }
                else if (div_lead[index] >= k_lag - overlap && div_lead[index] <= k_lag)
                {
                    if (!temp_indiv.empty()) // block starts with div value < in the desired range
                    {
                        for (auto x : temp_indiv)
                        {
                            block_info[x] = block_id;
                            group_info[x] = group_id; // odd id assigned
                        }
                        group_id += 1;
                        block_info[ppa_lead[index]] = block_id;
                        group_info[ppa_lead[index]] = group_id; // even id assigned
                        start = false;
                        temp_indiv.clear();
                    }
                    else // block starts with div value in the desired range
                    {
                        group_id += 1;
                        block_info[ppa_lead[index]] = block_id;
                        group_info[ppa_lead[index]] = group_id; // even id assigned
                        start = false;
                    }
                }
                else // No block
                {
                    temp_indiv.clear();
                }
            }
            else
            {
                if (div_lead[index] < k_lag - overlap)
                {
                    group_id += 1;
                    block_info[ppa_lead[index]] = block_id;
                    group_info[ppa_lead[index]] = group_id;
                }
                else if (div_lead[index] >= k_lag - overlap && div_lead[index] <= k_lag)
                {
                    if (group_id % 2 != 0)
                        group_id += 1;

                    block_info[ppa_lead[index]] = block_id;
                    group_info[ppa_lead[index]] = group_id;
                }
                else
                {
                    block_info[ppa_lead[index]] = block_id;
                    group_info[ppa_lead[index]] = group_id;
                    block_id += 1;
                    group_id = 1;
                    start = true;
                }
            }
        }
        return;
    }

    for (auto i = 0; i < M; ++i)
    {
        auto div_index = M - 1 - i;
        if (div_index >= 0)
        {
            if (start)
            {
                // block boundary possible start
                if (div_lead[div_index] < k_lag - overlap)
                {
                    temp_indiv.push_back(ppa_lead[div_index]);
                }
                else if (div_lead[div_index] >= k_lag - overlap && div_lead[div_index] <= k_lag)
                {
                    if (!temp_indiv.empty()) // block starts with div value < in the desired range
                    {
                        for (auto x : temp_indiv)
                        {
                            block_info[x] = block_id;
                            group_info[x] = group_id; // odd id assigned
                        }
                        group_id += 1;
                        block_info[ppa_lead[div_index]] = block_id;
                        group_info[ppa_lead[div_index]] = group_id; // even id assigned
                        start = false;
                        temp_indiv.clear();
                    }
                    else // block starts with div value in the desired range
                    {
                        group_id += 1;
                        block_info[ppa_lead[div_index]] = block_id;
                        group_info[ppa_lead[div_index]] = group_id; // even id assigned
                        start = false;
                    }
                }
                else // No block
                {
                    temp_indiv.clear();
                }
            }
            else
            {
                if (div_lead[div_index] < k_lag - overlap)
                {
                    if (group_id % 2 == 0) // odd id assigned
                        group_id += 1;
                    block_info[ppa_lead[div_index]] = block_id;
                    group_info[ppa_lead[div_index]] = group_id;
                }
                else if (div_lead[div_index] >= k_lag - overlap && div_lead[div_index] <= k_lag)
                {
                    if (group_id % 2 != 0) // even id assigned
                        group_id += 1;

                    block_info[ppa_lead[div_index]] = block_id;
                    group_info[ppa_lead[div_index]] = group_id;
                }
                else
                {
                    if (group_id % 2 == 0)
                        group_id += 1;

                    block_info[ppa_lead[div_index]] = block_id;
                    group_info[ppa_lead[div_index]] = group_id;
                    block_id += 1;
                    group_id = 1;
                    start = true;
                }
            }
        }

        // update ppa and divergence arrays
        auto index = ppa_lead[i];
        auto match_start = div_lead[i];
        auto allele = X_lead[index] - '0';

        if (match_start > p)
            p = match_start;
        if (match_start > q)
            q = match_start;

        if (allele == 0)
        {
            a.push_back(index);
            d.push_back(p);
            p = 0;
        }
        else
        {
            b.push_back(index);
            e.push_back(q);
            q = 0;
        }
    }

    if (k_lead < N)
    {
        ppa_lead = a;
        ppa_lead.insert(ppa_lead.end(), b.begin(), b.end());

        div_lead = d;
        div_lead.insert(div_lead.end(), e.begin(), e.end());
    }
    // End of Leading PBWT!
}

void laggingAlgo3(std::vector<char> &X_lag, double L, std::vector<int32_t> &ppa_lag,
                  std::vector<int32_t> &div_lag, int32_t k_lag, int32_t N,
                  std::vector<int32_t> &block_info, std::vector<int32_t> &group_info,
                  std::vector<std::string> &individuals, std::vector<std::string> &site,
                  std::map<int32_t, double> &genDist, int32_t k_lead, int32_t overlap,
                  std::vector<std::int32_t> &div_lead, std::vector<std::int32_t> &ppa_lead)
{
    std::vector<int32_t> a, b, d, e, ma, mb;
    auto q = k_lag + 1;
    auto p = q;

    for (auto i = 0; i < ppa_lag.size(); i++)
    {

        // Find matches ending in a site
        if (genDist[div_lag[i]] > genDist[k_lag] - L)
        {
            if (!ma.empty() && !mb.empty())
            {
                if (overlap == 0)
                {
                    findAltMatchesLessInfo(ma, mb, block_info, group_info, k_lag, individuals, site, k_lead);
                }
                else
                {
                    findAltMatchesLessInfo(ma, mb, block_info, group_info, k_lag, individuals, site, k_lead,
                                           overlap, div_lead, ppa_lead);
                }
            }
            ma.clear();
            mb.clear();
        }

        // update ppa and divergence arrays
        auto index = ppa_lag[i];
        auto match_start = div_lag[i];
        auto allele = X_lag[index] - '0';

        if (match_start > p)
            p = match_start;
        if (match_start > q)
            q = match_start;

        if (allele == 0)
        {
            a.push_back(index);
            d.push_back(p);
            p = 0;
            ma.push_back(index);
        }
        else
        {
            b.push_back(index);
            e.push_back(q);
            q = 0;
            mb.push_back(index);
        }

    } // iterating over ppa-values

    if (!ma.empty() && !mb.empty())
    {

        if (overlap == 0)
        {
            findAltMatchesLessInfo(ma, mb, block_info, group_info, k_lag, individuals, site, k_lead);
        }
        else
        {
            findAltMatchesLessInfo(ma, mb, block_info, group_info, k_lag, individuals, site, k_lead,
                                   overlap, div_lead, ppa_lead);
        }
    }

    if (k_lag < N)
    {
        ppa_lag = a;
        ppa_lag.insert(ppa_lag.end(), b.begin(), b.end());

        div_lag = d;
        div_lag.insert(div_lag.end(), e.begin(), e.end());
    }
    // std::cout << "End of Lagginb PBWT!\n";
}

int32_t findStartPosMatch(int32_t indv1, int32_t indv2, const std::vector<int32_t> &div, const std::vector<int32_t> &ppa)
{
    auto it1 = find(ppa.begin(), ppa.end(), indv1);
    auto it2 = find(ppa.begin(), ppa.end(), indv2);

    int32_t indv1_index = distance(ppa.begin(), it1);
    int32_t indv2_index = distance(ppa.begin(), it2);

    int32_t start = 0, end = 0;
    if (indv1_index < indv2_index)
    {
        start = indv1_index;
        end = indv2_index;
    }
    else
    {
        start = indv2_index;
        end = indv1_index;
    }

    int32_t start_pos = *max_element(div.begin() + start + 1, div.begin() + end + 1);
    return start_pos;
}

// Report Alternating Matches for Strict Boundary
void findAltMatchesLessInfo(std::vector<int32_t> &ma, std::vector<int32_t> &mb,
                            std::vector<int32_t> &block_info, std::vector<int32_t> &group_info, int32_t k_lag,
                            std::vector<std::string> &individuals, std::vector<std::string> &site,
                            int32_t k_lead)
{

    int32_t target_id1, target_id2, id1, id2;
    std::string actual_id1, actual_id2;
    for (auto i = 0; i < ma.size(); ++i)
    {
        id1 = ma[i];
        actual_id1 = individuals[id1].substr(0, individuals[id1].size() - 1);

        if (id1 % 2 == 0)
            target_id1 = id1 + 1;
        else
            target_id1 = id1 - 1;

        for (auto j = 0; j < mb.size(); ++j)
        {
            id2 = mb[j];
            actual_id2 = individuals[id2].substr(0, individuals[id2].size() - 1);

            // same individual different haplotype
            if (actual_id1 == actual_id2)
                continue;

            if (id2 % 2 == 0)
                target_id2 = id2 + 1;
            else
                target_id2 = id2 - 1;

            if (block_info[id1] == block_info[target_id2] && block_info[id1] > 0)
            {
                if (group_info[id1] * group_info[target_id2] < 0)
                {
                    // indv2 switched
                    std::cout << k_lag - 1 << ", " << k_lag << ", " << site[k_lag - 1] << ", " << site[k_lag]
                              << ", " << actual_id1 << ", " << actual_id2 << ", "
                              << individuals[id1].back() << ", " << individuals[id2].back() << ", "
                              << individuals[id1].back() << ", " << individuals[target_id2].back() << ", A\n";
                }
            }
            if (block_info[target_id1] == block_info[id2] && block_info[id2] > 0)
            {
                if (group_info[target_id1] * group_info[id2] < 0)
                {
                    // indv1 switched
                    std::cout << k_lag - 1 << ", " << k_lag << ", " << site[k_lag - 1] << ", " << site[k_lag]
                              << ", " << actual_id2 << ", " << actual_id1 << ", "
                              << individuals[id2].back() << ", " << individuals[id1].back() << ", "
                              << individuals[id2].back() << ", " << individuals[target_id1].back() << ", A\n";
                }
            }
            if (block_info[target_id1] == block_info[target_id2] && block_info[target_id2] > 0)
            {
                if (group_info[target_id1] * group_info[target_id2] < 0)
                {
                    // Double switched
                    std::cout << k_lag - 1 << ", " << k_lag << ", " << site[k_lag - 1]
                              << ", " << site[k_lag] << ", " << actual_id1 << ", " << actual_id2 << ", "
                              << individuals[id1].back() << ", " << individuals[id2].back() << ", "
                              << individuals[target_id1].back() << ", " << individuals[target_id2].back() << ", D\n";
                }
            }
        }
    }
}

// Report Alternating Matches for Overlap
void findAltMatchesLessInfo(std::vector<int32_t> &ma, std::vector<int32_t> &mb,
                            std::vector<int32_t> &block_info, std::vector<int32_t> &group_info, int32_t k_lag,
                            std::vector<std::string> &individuals, std::vector<std::string> &site,
                            int32_t k_lead, int32_t overlap, std::vector<int32_t> &div_lead,
                            std::vector<int32_t> &ppa_lead)
{

    int32_t target_id1, target_id2, id1, id2;
    std::string actual_id1, actual_id2;
    for (auto i = 0; i < ma.size(); ++i)
    {
        id1 = ma[i];
        actual_id1 = individuals[id1].substr(0, individuals[id1].size() - 1);

        if (id1 % 2 == 0) // find complement haplotype
            target_id1 = id1 + 1;
        else
            target_id1 = id1 - 1;

        for (auto j = 0; j < mb.size(); ++j)
        {
            id2 = mb[j];
            actual_id2 = individuals[id2].substr(0, individuals[id2].size() - 1);

            // same individual different haplotype
            if (actual_id1 == actual_id2)
                continue;

            if (id2 % 2 == 0) // find complement haplotype
                target_id2 = id2 + 1;
            else
                target_id2 = id2 - 1;

            int32_t leading_start_pos = 0;

            // Indv2 switched
            if (block_info[id1] == block_info[target_id2] && block_info[id1] > 0)
            {
                bool match_exists = true;
                if (group_info[id1] != group_info[target_id2])
                {
                    if (group_info[id1] % 2 == 0 && group_info[target_id2] % 2 != 0 && (group_info[id1] == group_info[target_id2] + 1)) // odd --> even, probable case that match doesn't exist
                    {
                        leading_start_pos = findStartPosMatch(id1, target_id2, div_lead, ppa_lead);
                        if (leading_start_pos < k_lag - overlap || leading_start_pos > k_lag)
                        {
                            match_exists = false;
                        }
                    }
                    else if (group_info[target_id2] % 2 == 0 && group_info[id1] % 2 != 0 && (group_info[target_id2] == group_info[id1] + 1)) // odd --> even, probable case that match doesn't exist
                    {
                        leading_start_pos = findStartPosMatch(id1, target_id2, div_lead, ppa_lead);
                        if (leading_start_pos < k_lag - overlap || leading_start_pos > k_lag)
                        {
                            match_exists = false;
                        }
                    }
                }
                else if (group_info[id1] == group_info[target_id2] && group_info[id1] % 2 != 0 && group_info[target_id2] % 2 != 0) // belong to same odd group, match doesn't exist
                {
                    match_exists = false;
                }
                else
                {
                    match_exists = true;
                }

                if (match_exists)
                {
                    std::cout << k_lag - 1 << "," << k_lag << ", " << site[k_lag - 1] << ", " << site[k_lag]
                              << ", " << actual_id1 << ", " << actual_id2 << ", "
                              << individuals[id1].back() << ", " << individuals[id2].back() << ", "
                              << individuals[id1].back() << ", " << individuals[target_id2].back() << ", A\n";
                }
            }

            // Indv1 switched
            if (block_info[target_id1] == block_info[id2] && block_info[id2] > 0)
            {
                bool match_exists = true;

                if (group_info[target_id1] != group_info[id2])
                {
                    if (group_info[target_id1] % 2 == 0 && group_info[id2] % 2 != 0 && (group_info[target_id1] == group_info[id2] + 1)) // odd --> even, probable case that match doesn't exist
                    {
                        leading_start_pos = findStartPosMatch(target_id1, id2, div_lead, ppa_lead);
                        if (leading_start_pos < k_lag - overlap || leading_start_pos > k_lag)
                        {
                            match_exists = false;
                        }
                    }
                    else if (group_info[id2] % 2 == 0 && group_info[target_id1] % 2 != 0 && (group_info[id2] == group_info[target_id1] + 1)) // odd --> even, probable case that match doesn't exist exist
                    {
                        leading_start_pos = findStartPosMatch(target_id1, id2, div_lead, ppa_lead);
                        if (leading_start_pos < k_lag - overlap || leading_start_pos > k_lag)
                        {
                            match_exists = false;
                        }
                    }
                }
                else if (group_info[target_id1] == group_info[id2] && group_info[target_id1] % 2 != 0 && group_info[id2] % 2 != 0) // belong to same odd group, match doesn't exist
                {
                    match_exists = false;
                }
                else
                {
                    match_exists = true;
                }

                if (match_exists)
                {
                    std::cout << k_lag - 1 << "," << k_lag << ", " << site[k_lag - 1] << ", " << site[k_lag]
                              << ", " << actual_id2 << ", " << actual_id1 << ", "
                              << individuals[id2].back() << ", " << individuals[id1].back() << ", "
                              << individuals[id2].back() << ", " << individuals[target_id1].back() << ", A\n";
                }
            }

            // // Double Switched
            // //     Uncomment the block of code below to output Double Switches
            if (block_info[target_id1] == block_info[target_id2] && block_info[target_id2] > 0)
            {
                bool match_exists = true;

                if (group_info[target_id1] != group_info[target_id2])
                {
                    if (group_info[target_id1] % 2 == 0 && group_info[target_id2] % 2 != 0 &&
                        (group_info[target_id1] == group_info[target_id2] + 1)) // odd followed by even by just one more, probable case match doesn't exist
                    {
                        leading_start_pos = findStartPosMatch(target_id1, target_id2, div_lead, ppa_lead);
                        if (leading_start_pos < k_lag - overlap || leading_start_pos > k_lag)
                        {
                            match_exists = false;
                        }
                    }
                    else if (group_info[target_id2] % 2 == 0 && group_info[target_id1] % 2 != 0 &&
                             (group_info[target_id2] == group_info[target_id1] + 1)) // odd followed by even by just one more, probable case match doesn't exist
                    {
                        leading_start_pos = findStartPosMatch(target_id1, target_id2, div_lead, ppa_lead);
                        if (leading_start_pos < k_lag - overlap || leading_start_pos > k_lag)
                        {
                            match_exists = false;
                        }
                    }
                }
                else if (group_info[target_id1] == group_info[target_id2] &&
                         group_info[target_id1] % 2 != 0 && group_info[target_id2] % 2 != 0)
                {
                    match_exists = false;
                }
                else
                {
                    match_exists = true;
                }

                if (match_exists)
                {
                    std::cout << k_lag - 1 << "," << k_lag << ", " << site[k_lag - 1] << ", "
                              << site[k_lag] << ", " << actual_id1 << ", " << actual_id2 << ", "
                              << individuals[id1].back() << ", " << individuals[id2].back() << ", "
                              << individuals[target_id1].back() << ", " << individuals[target_id2].back() << ", D\n";
                }
            }
        }
    }
}

// For strict Boundary
void updateBlockGroup(std::vector<int32_t> &div_lead, std::vector<int32_t> &ppa_lead,
                      int32_t N, std::vector<int32_t> &local_block_info,
                      std::vector<int32_t> &local_group_info, std::vector<std::string> &individuals, std::vector<std::string> &site, int32_t k_lag, int32_t k_lead)
{

    int32_t M = ppa_lead.size();
    bool find_block_start = true;
    std::vector<int32_t> curr_indiv;
    int32_t block_id = 1;

    for (auto i = 0; i < M; ++i)
    {
        auto div_index = M - 1 - i;
        if (div_index >= 0)
        {
            if (find_block_start)
            {
                // block boundary possible start
                if (div_lead[div_index] < k_lag)
                {
                    curr_indiv.push_back(ppa_lead[div_index]);
                }
                else if (div_lead[div_index] == k_lag)
                {
                    // VALID BLOCK
                    curr_indiv.push_back(ppa_lead[div_index]);
                    for (auto indiv : curr_indiv)
                    {
                        local_block_info[indiv] = block_id;
                        local_group_info[indiv] = 1;
                    }

                    find_block_start = false;
                    curr_indiv.clear();
                }
                else
                {
                    curr_indiv.clear();
                    // continue;
                }
            }
            else
            {
                // block boundary continues
                if (div_lead[div_index] < k_lag)
                    curr_indiv.push_back(ppa_lead[div_index]);
                else
                {
                    // block boundary ends
                    curr_indiv.push_back(ppa_lead[div_index]);
                    for (auto indiv : curr_indiv)
                    {
                        local_block_info[indiv] = block_id;
                        local_group_info[indiv] = -1;
                    }
                    find_block_start = true;
                    curr_indiv.clear();
                    block_id += 1;
                }
            }
        }
    }
}

// For overlap
void updateBlockGroup(std::vector<int32_t> &div_lead, std::vector<int32_t> &ppa_lead,
                      int32_t N, std::vector<int32_t> &local_block_info, std::vector<int32_t> &local_group_info,
                      std::vector<std::string> &individuals, std::vector<std::string> &site,
                      int32_t k_lag, int32_t k_lead, int32_t overlap)
{

    int32_t M = ppa_lead.size();
    std::vector<int32_t> temp_indiv;
    int32_t block_id = 1, group_id = 1;
    bool start = true;

    for (auto div_index = M - 1; div_index >= 0; --div_index)
    {
        if (start)
        {
            // block boundary possible start
            if (div_lead[div_index] < k_lag - overlap)
            {
                temp_indiv.push_back(ppa_lead[div_index]);
            }
            else if (div_lead[div_index] >= k_lag - overlap && div_lead[div_index] <= k_lag)
            {
                if (!temp_indiv.empty()) // block starts with div value < in the desired range
                {
                    for (auto x : temp_indiv)
                    {
                        local_block_info[x] = block_id;
                        local_group_info[x] = group_id; // odd id assigned
                    }
                    group_id += 1;
                    local_block_info[ppa_lead[div_index]] = block_id;
                    local_group_info[ppa_lead[div_index]] = group_id; // even id assigned
                    start = false;
                    temp_indiv.clear();
                }
                else // block starts with div value in the desired range
                {
                    group_id += 1;
                    local_block_info[ppa_lead[div_index]] = block_id;
                    local_group_info[ppa_lead[div_index]] = group_id; // even id assigned
                    start = false;
                }
            }
            else // No block
            {
                temp_indiv.clear();
            }
        }
        else
        {
            if (div_lead[div_index] < k_lag - overlap)
            {
                if (group_id % 2 == 0)
                    group_id += 1;
                local_block_info[ppa_lead[div_index]] = block_id;
                local_group_info[ppa_lead[div_index]] = group_id;
            }
            else if (div_lead[div_index] >= k_lag - overlap && div_lead[div_index] <= k_lag)
            {
                if (group_id % 2 != 0)
                    group_id += 1;

                local_block_info[ppa_lead[div_index]] = block_id;
                local_group_info[ppa_lead[div_index]] = group_id;
            }
            else
            {
                if (group_id % 2 == 0)
                    group_id += 1;
                local_block_info[ppa_lead[div_index]] = block_id;
                local_group_info[ppa_lead[div_index]] = group_id;
                block_id += 1;
                group_id = 1;
                start = true;
            }
        }
    }
}

std::map<int32_t, double> getGenDist(std::string genMapPath)
{
    std::cout << "** Loading Genetic map...\n";
    std::map<int32_t, double> genDist;
    std::string line = "";
    std::ifstream inFile(genMapPath);

    int32_t site = 0;
    while (getline(inFile, line))
    {
        std::istringstream is(line);
        std::string word;
        int32_t pos = 0;
        // int32_t key = 0;
        double value = 0.0;
        while (is >> word)
        {
            if (pos == 2)
                value = std::stod(word);
            pos += 1;
        }
        genDist[site] = value;
        site += 1;
    }
    std::cout << "** sucessfully loaded genetic map!\n";
    return genDist;
}
