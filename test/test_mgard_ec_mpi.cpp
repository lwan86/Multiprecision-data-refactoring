#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <queue>
#include <unordered_map>
#include <list>
#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>

#include "../include/Decomposer/Decomposer.hpp"
#include "../include/Interleaver/Interleaver.hpp"
#include "../include/BitplaneEncoder/BitplaneEncoder.hpp"
#include "../include/ErrorEstimator/ErrorEstimator.hpp"
#include "../include/ErrorCollector/ErrorCollector.hpp"
#include "../include/LosslessCompressor/LevelCompressor.hpp"
#include "../include/RefactorUtils.hpp"

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

using namespace ROCKSDB_NAMESPACE;

#include <erasurecode.h>
#include <erasurecode_helpers.h>
#include <config_liberasurecode.h>
#include <erasurecode_stdinc.h>
#include <erasurecode_version.h>


#include <mpi.h>

#include <adios2.h>

struct UnitErrorGain{
    double unit_error_gain;
    int level;
    UnitErrorGain(double u, int l) : unit_error_gain(u), level(l) {}
};
struct CompareUniteErrorGain{
    bool operator()(const UnitErrorGain& u1, const UnitErrorGain& u2){
        return u1.unit_error_gain < u2.unit_error_gain;
    }
};

template <class T>
std::vector<std::tuple<uint64_t, uint64_t>> calculate_retrieve_order(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, double tolerance, std::vector<uint8_t>& index, MDR::MaxErrorEstimatorOB<T> error_estimator) 
{
    size_t num_levels = level_sizes.size();
    std::vector<std::tuple<uint64_t, uint64_t>> retrieve_order;
    double accumulated_error = 0;
    for(size_t i = 0; i < num_levels; i++)
    {
        accumulated_error += error_estimator.estimate_error(level_errors[i][index[i]], i);
    }
    std::priority_queue<UnitErrorGain, std::vector<UnitErrorGain>, CompareUniteErrorGain> heap;
    // identify minimal level
    double min_error = accumulated_error;
    for(size_t i = 0; i < num_levels; i++)
    {
        min_error -= error_estimator.estimate_error(level_errors[i][index[i]], i);
        min_error += error_estimator.estimate_error(level_errors[i].back(), i);
        // fetch the first component if index is 0
        if(index[i] == 0){
            //retrieve_sizes[i] += level_sizes[i][index[i]];
            accumulated_error -= error_estimator.estimate_error(level_errors[i][index[i]], i);
            accumulated_error += error_estimator.estimate_error(level_errors[i][index[i] + 1], i);
            retrieve_order.push_back(std::make_tuple(i, static_cast<uint64_t>(index[i])));
            index[i] ++;
            //std::cout << i;
        }
        // push the next one
        if(index[i] != level_sizes[i].size())
        {
            double error_gain = error_estimator.estimate_error_gain(accumulated_error, level_errors[i][index[i]], level_errors[i][index[i]+1], i);
            heap.push(UnitErrorGain(error_gain/level_sizes[i][index[i]], i));
        }
        //std::cout << i;
        //std::cout << i << ", " << retrieve_sizes[i] << std::endl;
        if(min_error < tolerance)
        {
            // the min error of first 0~i levels meets the tolerance
            num_levels = i+1;
            break;
        }
    }

    bool tolerance_met = accumulated_error < tolerance;
    while((!tolerance_met) && (!heap.empty()))
    {
        auto unit_error_gain = heap.top();
        heap.pop();
        int i = unit_error_gain.level;
        int j = index[i];
        //retrieve_sizes[i] += level_sizes[i][j];
        accumulated_error -= error_estimator.estimate_error(level_errors[i][j], i);
        accumulated_error += error_estimator.estimate_error(level_errors[i][j+1], i);
        if(accumulated_error < tolerance)
        {
            tolerance_met = true;
        }
        //std::cout << i << ", " << +index[i] << std::endl;
        retrieve_order.push_back(std::make_tuple(i, static_cast<uint64_t>(index[i])));
        index[i]++;
        if(index[i] != level_sizes[i].size())
        {
            double error_gain = error_estimator.estimate_error_gain(accumulated_error, level_errors[i][index[i]], level_errors[i][index[i]+1], i);
            heap.push(UnitErrorGain(error_gain/level_sizes[i][index[i]], i));
        }
        //std::cout << i;
        //std::cout << i << ", " << retrieve_sizes[i] << std::endl;
    }
    //std::cout << std::endl;
    //std::cout << "Requested tolerance = " << tolerance << ", estimated error = " << accumulated_error << std::endl;
    return retrieve_order;
}

void decomposition(const size_t ndim, const size_t rank, const int *decomp,
                    int *pos, bool rowmajor)
{
    if (rowmajor)
    {
        size_t prod = 1;
        for (int k = static_cast<int>(ndim) - 1; k >= 0; --k)
        {
            if (decomp[k] == 0)
            {
                pos[k] = -1;
                continue;
            }
            pos[k] = rank / prod % decomp[k];
            prod *= decomp[k];
        }
    }
    else
    {
        size_t prod = 1;
        for (size_t k = 0; k < ndim; ++k)
        {
            if (decomp[k] == 0)
            {
                pos[k] = -1;
                continue;
            }
            pos[k] = rank / prod % decomp[k];
            prod *= decomp[k];
        }
    }

}

int main(int argc, char *argv[])
{
    std::string inputFileName;
    size_t dataTiers = 0;
    std::vector<double> dataTiersTolerance;

    std::vector<int> dataTiersECParam_k;
    std::vector<int> dataTiersECParam_m;
    std::vector<int> dataTiersECParam_w;
    std::string ECBackendName = "null";
    size_t total_mgard_levels = 0;
    size_t num_bitplanes = 0;

    size_t procGroups = 1;
    size_t iterSteps;

    ec_backend_id_t backendID;

    size_t totalDimensions;
    std::vector<int> decompositions;
    std::list<size_t> fixedDimensions; 

    for (size_t i = 0; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "-i" || arg == "--input")
        {
            if (i+1 < argc)
            {
                inputFileName = argv[i+1];
            }
            else
            {
                std::cerr << "--input option requires one argument." << std::endl;
                return 1;
            }            
        }  
        else if (arg == "-t" || arg == "--tiers")
        {
            if (i+1 < argc)
            {
                dataTiers = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--tiers option requires [# of tiers] to be set first." << std::endl;
                return 1;
            }            
            if (dataTiers)
            {
                if (i+1+dataTiers < argc)
                { 
                    for (size_t j = i+2; j < i+2+dataTiers; j++)
                    {
                        dataTiersTolerance.push_back(atof(argv[j]));
                    } 
                    for (size_t j = i+2+dataTiers; j < i+2+dataTiers*2; j++)
                    {
                        dataTiersECParam_k.push_back(atoi(argv[j]));
                    }
                    for (size_t j = i+2+dataTiers*2; j < i+2+dataTiers*3; j++)
                    {
                        dataTiersECParam_m.push_back(atoi(argv[j]));
                    }
                    for (size_t j = i+2+dataTiers*3; j < i+2+dataTiers*4; j++)
                    {
                        dataTiersECParam_w.push_back(atoi(argv[j]));
                    }
                }
                else
                {
                    std::cerr << "--tiers option requires [# of tiers]*4 arguments." << std::endl;
                    return 1;
                } 
            }
            else
            {
                std::cerr << "--[# of tiers] needs to be greater than zero." << std::endl;
                return 1;
            }
        }
        else if (arg == "-b" || arg == "--backend")
        {
            if (i+1 < argc)
            {
                ECBackendName = argv[i+1];
            }
            else
            {
                std::cerr << "--backend option requires one argument." << std::endl;
                return 1;
            }            
        } 
        else if (arg == "-l" || arg == "--levels")
        {
            if (i+1 < argc)
            {
                total_mgard_levels = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--levels option requires one argument." << std::endl;
                return 1;
            }            
        }  
        else if (arg == "-p" || arg == "--planes")
        {
            if (i+1 < argc)
            {
                num_bitplanes = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--planes option requires one argument." << std::endl;
                return 1;
            }            
        }
        else if (arg == "-g" || arg == "--groups")
        {
            if (i+1 < argc)
            {
                procGroups = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--groups option requires one argument." << std::endl;
                return 1;
            }            
        }
        else if (arg == "-s" || arg == "--steps")
        {
            if (i+1 < argc)
            {
                iterSteps = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--steps option requires one argument." << std::endl;
                return 1;
            }            
        }
        else if (arg == "-d" || arg == "--decomposition")
        {
            if (i+1 < argc)
            {
                totalDimensions = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--decomposition option requires [# of dimensions] to be set first." << std::endl;
                return 1;
            }            
            if (totalDimensions)
            {
                if (i+1+totalDimensions < argc)
                {   
                    bool hasFixedDimension = false;
                    size_t fixedDimensionCount = 0;
                    for (size_t j = i+2; j < i+2+totalDimensions; j++)
                    {
                        int ranksPerDimension = atoi(argv[j]);
                        if (ranksPerDimension == 0)
                        {
                            if (!hasFixedDimension)
                            {
                                hasFixedDimension = true;
                            }
                            fixedDimensionCount++;
                        }
                        decompositions.push_back(ranksPerDimension);
                    }
                    if (hasFixedDimension)
                    {
                        if (i+1+totalDimensions+fixedDimensionCount < argc)
                        {
                            for (size_t j = i+2+totalDimensions; j < i+2+totalDimensions+fixedDimensionCount; j++)
                            {
                                size_t fixedDim = atoi(argv[j]);
                                fixedDimensions.push_back(fixedDim);
                            }
                        }
                        
                    }     
                }
                else
                {
                    std::cerr << "--decomposition option requires [# of dimensions] arguments." << std::endl;
                    return 1;
                } 
            }
            else
            {
                std::cerr << "--[# of dimensions] needs to be greater than zero." << std::endl;
                return 1;
            }
        }
    }
    int world_rank, world_size, rank, size;
    MPI_Init(&argc, &argv);
    
    MPI_Comm mpiComm;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int mpiGroup = world_rank%procGroups;
    MPI_Comm_split(MPI_COMM_WORLD, mpiGroup, world_rank, &mpiComm);
    MPI_Comm_rank(mpiComm, &rank);
    MPI_Comm_size(mpiComm, &size);
    std::cout << "global rank " << world_rank << ": I am rank " << rank << " of group " << mpiGroup << ", which has " << size << " ranks in total" << std::endl;
    
    if (ECBackendName == "flat_xor_hd")
    {
        backendID = EC_BACKEND_FLAT_XOR_HD;
    }
    else if (ECBackendName == "jerasure_rs_vand")
    {
        backendID = EC_BACKEND_JERASURE_RS_VAND;
    }
    else if (ECBackendName == "jerasure_rs_cauchy")
    {
        backendID = EC_BACKEND_JERASURE_RS_CAUCHY;
    }
    else if (ECBackendName == "isa_l_rs_vand")
    {
        backendID = EC_BACKEND_ISA_L_RS_VAND;
    }
    else if (ECBackendName == "isa_l_rs_cauchy")
    {
        backendID = EC_BACKEND_ISA_L_RS_CAUCHY;
    }
    else if (ECBackendName == "shss")
    {
        backendID = EC_BACKEND_SHSS;
    }
    else if (ECBackendName == "liberasurecode_rs_vand")
    {
        backendID = EC_BACKEND_LIBERASURECODE_RS_VAND;
    }
    else if (ECBackendName == "libphazr")
    {
        backendID = EC_BACKEND_LIBPHAZR;
    }
    else if (ECBackendName == "null")
    {
        backendID = EC_BACKEND_NULL;
    }
    else 
    {
        std::cerr << "the specified EC backend is not supported!" << std::endl;
        return 1;
    }
    //std::cout << "hello!" << std::endl;

    DB* db;
    Options options;

    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO reader_io = adios.DeclareIO("ReaderIO");
    adios2::Engine reader_engine =
           reader_io.Open(inputFileName, adios2::Mode::Read);
    const std::map<std::string, adios2::Params> allVariables =
           reader_io.AvailableVariables();
    
    std::unordered_map<std::string, std::string> allVarTypeMap;
    std::unordered_map<std::string, std::vector<std::vector<float>>> allVarDataMap;
    std::unordered_map<std::string, std::vector<std::vector<size_t>>> allVarStartMap;
    std::unordered_map<std::string, std::vector<std::vector<size_t>>> allVarCountMap;
    
    for (const auto variablePair : allVariables)
    {
        std::string variableName;
        std::string variableType;
        variableName = variablePair.first;
        variableType = variablePair.second.at("Type");  
        std::vector<int> variableDecomposition;
        
        if (variableType == "float")
        {
            auto variable = reader_io.InquireVariable<float>(variableName);
            size_t spaceDimensions = variable.Shape().size();
            if (totalDimensions != spaceDimensions)
            {
                continue;
            }
            if (world_rank == 0)
            {
                std::cout << "read " << spaceDimensions << "D variable " << variableName << std::endl;
            }
            if (spaceDimensions == 1)
            {
                variableDecomposition.push_back(world_size);
            }
            else
            {
                variableDecomposition = decompositions;
            }
            std::vector<int> pos(spaceDimensions);
            decomposition(spaceDimensions, world_rank, variableDecomposition.data(), pos.data(), true);
            
            std::vector<size_t> variablePerRankStart;
            std::vector<size_t> variablePerRankCount;
            size_t variablePerRankSize = 1;
            for (size_t i = 0; i < spaceDimensions; i++)
            {
                if (variableDecomposition[i] == 0)
                {
                    // certain dimension picked fixed plane
                    if (pos[i] != -1)
                    {
                        std::cerr << "pos should be -1 on dimension " << i << "."<< std::endl;
                        return 1;
                    }
                    size_t count = 1;
                    size_t offset = fixedDimensions.front();
                    fixedDimensions.pop_front();
                    variablePerRankStart.push_back(offset);
                    variablePerRankCount.push_back(count);
                }
                else
                {
                    size_t count = variable.Shape()[i]/variableDecomposition[i];
                    size_t offset = count*pos[i];
                    if (pos[i] == variableDecomposition[i]-1 && pos[i] != 0)
                    {
                        count = variable.Shape()[i]-offset;
                    }
                    variablePerRankStart.push_back(offset);
                    variablePerRankCount.push_back(count);
                    variablePerRankSize *= count;
                }
            }
            variable.SetSelection({variablePerRankStart, variablePerRankCount});

            std::vector<float> variablePerRankData(variablePerRankSize);
            reader_engine.Get(variable, variablePerRankData.data(), adios2::Mode::Sync);
            allVarTypeMap[variableName] = variableType;
            allVarDataMap[variableName].push_back(variablePerRankData);
            allVarStartMap[variableName].push_back(variablePerRankStart);
            allVarCountMap[variableName].push_back(variablePerRankCount);
        }
        break;
    }
    reader_engine.Close();

 

    for (size_t s = 0; s < iterSteps; s++)
    {

        for (auto &item : allVarDataMap)
        {
            if (allVarTypeMap[item.first] == "float")
            {
                using T = float;
                using T_stream = uint32_t;

                std::string variableName = item.first;
                int blocksCount = 0;

                std::vector<std::vector<std::vector<uint64_t>>> allBlocksQueryTable;
                std::vector<std::vector<std::vector<uint8_t>>> allBlocksDataTiersValues;
                std::vector<std::vector<T>> allBlocksLevelErrorBounds;
                std::vector<std::vector<std::vector<double>>> allBlocksLevelSquaredErrors;
                std::vector<std::vector<uint8_t>> allBlocksStoppingIndices; 
                
                std::vector<uint64_t> allBlocksQueryTableContent;
                std::vector<T> allBlocksLevelErrorBoundsContent;
                std::vector<double> allBlocksLevelSquaredErrorsContent;
                std::vector<uint8_t> allBlocksStoppingIndicesContent;

                for (auto & block : item.second)
                {
                    size_t block_id = blocksCount;
                    auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
                    auto interleaver = MDR::DirectInterleaver<T>();
                    auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
                    auto compressor = MDR::AdaptiveLevelCompressor(32);

                    std::vector<T> level_error_bounds;
                    std::vector<std::vector<uint8_t*>> level_components;
                    std::vector<std::vector<uint32_t>> level_sizes;
                    std::vector<std::vector<double>> level_squared_errors;
                    std::vector<uint8_t> stopping_indices;
                    std::vector<T> variableData = block;
                    std::vector<uint32_t> dimensions(totalDimensions);

                    for (size_t i = 0; i < allVarCountMap[variableName][block_id].size(); i++)
                    {
                        dimensions[i] = allVarCountMap[variableName][block_id][i];
                    }
                    uint8_t target_level = total_mgard_levels-1;
                    uint8_t max_level = log2(*min_element(dimensions.begin(), dimensions.end())) - 1;
                    if(target_level > max_level)
                    {
                        std::cerr << "Target level is higher than " << max_level << std::endl;
                        return 1;
                    }
                    // decompose data hierarchically
                    decomposer.decompose(variableData.data(), dimensions, target_level);

                    auto level_dims = MDR::compute_level_dims(dimensions, target_level);
                    auto level_elements = MDR::compute_level_elements(level_dims, target_level);
                    std::vector<uint32_t> dims_dummy(dimensions.size(), 0);
                    MDR::SquaredErrorCollector<T> s_collector = MDR::SquaredErrorCollector<T>();
                    for(size_t i = 0; i <= target_level; i++)
                    {
                        const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
                        T * buffer = (T *) malloc(level_elements[i] * sizeof(T));
                        // extract level i component
                        interleaver.interleave(variableData.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));
                        // compute max coefficient as level error bound
                        T level_max_error = MDR::compute_max_abs_value(reinterpret_cast<T*>(buffer), level_elements[i]);
                        level_error_bounds.push_back(level_max_error);
                        // collect errors
                        // auto collected_error = s_collector.collect_level_error(buffer, level_elements[i], num_bitplanes, level_max_error);
                        // std::cout << collected_error.size() << std::endl;
                        // level_squared_errors.push_back(collected_error);
                        // encode level data
                        int level_exp = 0;
                        frexp(level_max_error, &level_exp);
                        std::vector<uint32_t> stream_sizes;
                        std::vector<double> level_sq_err;
                        auto streams = encoder.encode(buffer, level_elements[i], level_exp, num_bitplanes, stream_sizes, level_sq_err);
                        free(buffer);
                        level_squared_errors.push_back(level_sq_err);
                        // lossless compression
                        uint8_t stopping_index = compressor.compress_level(streams, stream_sizes);
                        stopping_indices.push_back(stopping_index);
                        // record encoded level data and size
                        level_components.push_back(streams);
                        level_sizes.push_back(stream_sizes);
                    }
                    std::vector<uint8_t> level_num_bitplanes(target_level+1, 0);
                    std::vector<std::vector<uint8_t>> dataTiersValues(dataTiers);
                    std::vector<std::vector<uint64_t>> queryTable;
                    for (size_t i = 0; i < dataTiersTolerance.size(); i++)
                    {
                        // if (world_rank == 0)
                        // {
                        //     std::cout << "tolerance: " << dataTiersTolerance[i] << std::endl;
                        // }
                        std::vector<uint8_t> oneDataTierValues;
                        uint64_t currentTierCopidedSize = 0;
                        
                        std::vector<std::vector<double>> level_abs_errors;
                        std::vector<std::vector<double>> level_errors;
                        target_level = level_error_bounds.size() - 1;
                        MDR::MaxErrorCollector<T> collector = MDR::MaxErrorCollector<T>();
                        for(size_t j = 0; j <= target_level; j++)
                        {
                            auto collected_error = collector.collect_level_error(NULL, 0, level_squared_errors[j].size(), level_error_bounds[j]);
                            level_abs_errors.push_back(collected_error);
                        }
                        level_errors = level_abs_errors;
                        auto prev_level_num_bitplanes(level_num_bitplanes);
                        auto estimator = MDR::MaxErrorEstimatorOB<T>(totalDimensions); 
                        auto retrieve_order = calculate_retrieve_order(level_sizes, level_errors, dataTiersTolerance[i], level_num_bitplanes, estimator);
                        for (size_t j = 0; j < retrieve_order.size(); j++)
                        {
                            uint64_t lid = std::get<0>(retrieve_order[j]);
                            uint64_t pid = std::get<1>(retrieve_order[j]);
                            // if (world_rank == 0)
                            // {
                            //     std::cout << lid << ", " << pid << std::endl;
                            // }
                            std::vector<uint8_t> onePieceValues(level_components[lid][pid], level_components[lid][pid]+level_sizes[lid][pid]);
                            oneDataTierValues.insert(oneDataTierValues.end(), onePieceValues.begin(), onePieceValues.end());
                            //uint64_t group = mpiGroup;
                            //uint64_t group_member_id = rank;
                            uint64_t tier_id = i;
                            uint64_t tier_size = level_sizes[lid][pid];
                            std::vector<uint64_t> row = {lid, pid, tier_id, currentTierCopidedSize, tier_size};
                            queryTable.push_back(row);
                            currentTierCopidedSize += level_sizes[lid][pid];
                        }
                        dataTiersValues.push_back(oneDataTierValues);
                    }
                    allBlocksDataTiersValues.push_back(dataTiersValues);
                    allBlocksQueryTable.push_back(queryTable);
                    allBlocksLevelErrorBounds.push_back(level_error_bounds);
                    allBlocksLevelSquaredErrors.push_back(level_squared_errors);
                    allBlocksStoppingIndices.push_back(stopping_indices);
                    blocksCount++;
                }

                //std::cout << "global rank "<< world_rank << " queryTable size: " << queryTable.size() << std::endl;

                
                for (size_t i = 0; i < allBlocksQueryTable.size(); i++)
                {
                    for (size_t j = 0; j < allBlocksQueryTable[i].size(); j++)
                    {
                        allBlocksQueryTableContent.insert(allBlocksQueryTableContent.end(), allBlocksQueryTable[i][j].begin(), allBlocksQueryTable[i][j].end());
                    }   
                }
                for (size_t i = 0; i < allBlocksLevelErrorBounds.size(); i++)
                {
                    allBlocksLevelErrorBoundsContent.insert(allBlocksLevelErrorBoundsContent.end(), allBlocksLevelErrorBounds[i].begin(), allBlocksLevelErrorBounds[i].end());
                }
                for (size_t i = 0; i < allBlocksLevelSquaredErrors.size(); i++)
                {
                    for (size_t j = 0; j < allBlocksLevelSquaredErrors[i].size(); j++)
                    {
                        allBlocksLevelSquaredErrorsContent.insert(allBlocksLevelSquaredErrorsContent.end(), allBlocksLevelSquaredErrors[i][j].begin(), allBlocksLevelSquaredErrors[i][j].end());
                    }
                }
                for (size_t i = 0; i < allBlocksStoppingIndices.size(); i++)
                {
                    allBlocksStoppingIndicesContent.insert(allBlocksStoppingIndicesContent.end(), allBlocksStoppingIndices[i].begin(), allBlocksStoppingIndices[i].end());
                }      

                int allBlocksQueryTableContentLength = allBlocksQueryTableContent.size();  
                int allBlocksLevelErrorBoundsContentLength = allBlocksLevelErrorBoundsContent.size();
                int allBlocksLevelSquaredErrorsContentLength = allBlocksLevelSquaredErrorsContent.size();
                int allBlocksStoppingIndicesContentLength = allBlocksStoppingIndicesContent.size();


                int *recvQueryTableContentLength = NULL;
                int *recvBlocksCount = NULL;
                if (world_rank == 0)
                {
                    recvQueryTableContentLength = (int *) malloc(world_size * sizeof(int));
                    recvBlocksCount = (int *) malloc(world_size * sizeof(int));
                }
                MPI_Gather(&allBlocksQueryTableContentLength, 1, MPI_INT, recvQueryTableContentLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&blocksCount, 1, MPI_INT, recvBlocksCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
                
                int *displs = NULL;
                uint64_t *totalRecvQueryTableContent = NULL;

                if (world_rank == 0)
                {
                    displs = (int *) malloc(world_size * sizeof(int));
                    int totalRecvQueryTableContentLength = 0;
                    for (size_t i = 0; i < world_size; i++)
                    {
                        totalRecvQueryTableContentLength += recvQueryTableContentLength[i];
                    }
                    //std::cout << totalRecvQueryTableContentLength << std::endl;
                    totalRecvQueryTableContent = (uint64_t *) malloc(totalRecvQueryTableContentLength);
                    displs[0] = 0;
                    for (size_t i = 1; i < world_size; i++)
                    {
                        displs[i] = displs[i-1]+recvQueryTableContentLength[i-1];
                    }
                    
                }
                MPI_Gatherv(allBlocksQueryTableContent.data(), allBlocksQueryTableContentLength, MPI_UINT64_T,
                            totalRecvQueryTableContent, recvQueryTableContentLength, displs, MPI_UINT64_T,
                            0, MPI_COMM_WORLD);      

                
                if (world_rank == 0)
                {
                    // Optimize RocksDB. This is the easiest way to get RocksDB to perform well
                    // options.IncreaseParallelism();
                    // options.OptimizeLevelStyleCompaction();
                    // // create the DB if it's not already present
                    // options.create_if_missing = true;
                    // open DB
                    //Status st = DB::Open(options, "meta.db", &db);
                    //assert(st.ok());
                    for (size_t i = 0; i < world_size; i++)
                    {
                        int oneRankBlocksCount = recvBlocksCount[i];
                        std::vector<uint64_t> oneRankQueryTableContent(totalRecvQueryTableContent+displs[i], totalRecvQueryTableContent+displs[i]+recvQueryTableContentLength[i]);
                        std::cout << "rank " << i << " has " << oneRankBlocksCount << " blocks!" << std::endl;
                        int oneRankQueryTableRows = oneRankQueryTableContent.size()/5;
                        int oneRankQueryTableRowsPerBlock = oneRankQueryTableRows/oneRankBlocksCount;
                        int rowID = 0;
                        int blockID = 0;
                        //std::vector<uint64_t> oneBlockQueryTableContent;
                        for (size_t j = 0; j < oneRankQueryTableContent.size(); j+=5)
                        {
                            std::cout << "block " << blockID << ": ";
                            for (size_t k = j; k < j+5; k++)
                            {
                                std::cout << oneRankQueryTableContent[k] << ", ";
                                //oneBlockQueryTableContent.push_back(oneRankQueryTableContent[k]);
                            }
                            std::cout << std::endl;
                            rowID++;
                            if (rowID%oneRankQueryTableRowsPerBlock == 0)
                            {
                                // std::ostringstream vts;
                                // std::copy(oneBlockQueryTableContent.begin(), oneBlockQueryTableContent.end()-1,
                                //         std::ostream_iterator<int>(vts, ", "));
                                // vts << oneBlockQueryTableContent.back();
                                // std::cout << vts.str() << std::endl;
                                //st = db->Put(WriteOptions(), variableName+":step:"+std::to_string(s)+":rank:"+std::to_string(i)+":block:"+std::to_string(blockID), vts.str());
                                //assert(st.ok());
                                //oneBlockQueryTableContent.clear();
                                blockID++;
                            }                            
                        }
                        

                    }
                    
                }
                          

                if (world_rank == 0)
                {
                    free(displs);
                    free(recvQueryTableContentLength);
                    free(recvBlocksCount);
                    free(totalRecvQueryTableContent);
                }
                
                
                
                // if (world_rank == 1)
                // {
                //     for (size_t i = 0; i < queryTable.size(); i++)
                //     {
                //         for (size_t j = 0; j < queryTable[i].size(); j++)
                //         {
                //             std::cout << queryTable[i][j] << ", ";
                //         }
                //         std::cout << std::endl;
                //     }
                // }
                
                
            }
            

            
        }

        
    }

    //delete db;
    MPI_Finalize();

}