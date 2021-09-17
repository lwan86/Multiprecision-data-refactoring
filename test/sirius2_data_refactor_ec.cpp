#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <queue>
#include <unordered_map>

#include "../include/Decomposer/Decomposer.hpp"
#include "../include/Interleaver/Interleaver.hpp"
#include "../include/BitplaneEncoder/BitplaneEncoder.hpp"
#include "../include/ErrorEstimator/ErrorEstimator.hpp"
#include "../include/ErrorCollector/ErrorCollector.hpp"
#include "../include/LosslessCompressor/LevelCompressor.hpp"
#include "../include/RefactorUtils.hpp"

#include <erasurecode.h>
#include <erasurecode_helpers.h>
#include <config_liberasurecode.h>
#include <erasurecode_stdinc.h>
#include <erasurecode_version.h>

//#include <liberasurecode/erasurecode_backend.h> 
//#include <liberasurecode/erasurecode_helpers_ext.h> 
//#include <liberasurecode/erasurecode_preprocessing.h> 

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
std::vector<std::tuple<uint32_t, uint32_t>> calculate_retrieve_order(const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, double tolerance, std::vector<uint8_t>& index, MDR::MaxErrorEstimatorOB<T> error_estimator) 
{
    size_t num_levels = level_sizes.size();
    std::vector<std::tuple<uint32_t, uint32_t>> retrieve_order;
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
            retrieve_order.push_back(std::make_tuple(i, static_cast<uint32_t>(index[i])));
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
        retrieve_order.push_back(std::make_tuple(i, static_cast<uint32_t>(index[i])));
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

int main(int argc, char *argv[])
{
    std::string inputFileName;
    size_t storageTiers = 0;
    std::vector<std::string> storageTiersPaths;
    std::vector<double> storageTiersPercentages;
    
    std::vector<int> storageTiersECParam_k;
    std::vector<int> storageTiersECParam_m;
    std::vector<int> storageTiersECParam_w;
    std::string ECBackendName = "null";
    size_t target_level = 0;
    size_t num_bitplanes = 0;

    ec_backend_id_t backendID;

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
                storageTiers = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--tiers option requires [# of tiers] to be set first." << std::endl;
                return 1;
            }            
            if (storageTiers)
            {
                if (i+1+storageTiers < argc)
                {
                    for (size_t j = i+2; j < i+2+storageTiers; j++)
                    {
                        char *abs_path;
                        abs_path = realpath(argv[j], NULL); 
                        std::string str(abs_path);
                        storageTiersPaths.push_back(str);
                    }  
                    for (size_t j = i+2+storageTiers; j < i+2+storageTiers*2; j++)
                    {
                        storageTiersPercentages.push_back(atof(argv[j]));
                    } 
                    for (size_t j = i+2+storageTiers*2; j < i+2+storageTiers*3; j++)
                    {
                        storageTiersECParam_k.push_back(atoi(argv[j]));
                    }
                    for (size_t j = i+2+storageTiers*3; j < i+2+storageTiers*4; j++)
                    {
                        storageTiersECParam_m.push_back(atoi(argv[j]));
                    }
                    for (size_t j = i+2+storageTiers*4; j < i+2+storageTiers*5; j++)
                    {
                        storageTiersECParam_w.push_back(atoi(argv[j]));
                    }
                }
                else
                {
                    std::cerr << "--tiers option requires [# of tiers]*5 arguments." << std::endl;
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
                target_level = atoi(argv[i+1]);
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
    } 

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

    adios2::ADIOS adios;
    adios2::IO writer_io = adios.DeclareIO("WriterIO");
    //std::vector<adios2::Engine> data_writer_engines(storageTiersPaths.size());
    std::unordered_map<std::string, adios2::Engine> data_writer_engines;
    std::unordered_map<std::string, adios2::Engine> parity_writer_engines;

    for (size_t i = 0; i < storageTiersPaths.size(); i++)
    {
        std::cout << "tier " << i << ": " << storageTiersPaths[i] << ", " << storageTiersPercentages[i] << ", " <<  storageTiersECParam_k[i] << ", " << storageTiersECParam_m[i] << ", " << storageTiersECParam_w[i] << std::endl;
        std::string fullDataPath = storageTiersPaths[i];
        //std::string refactoredDataFileName = "refactored.data.bp";
        if (!fullDataPath.empty() && fullDataPath.back() != '/')
        {
            fullDataPath += '/';
        }
        for (size_t j = 0; j < storageTiersECParam_k[i]; j++)
        {
            std::string idxStr = std::to_string(i)+"_"+std::to_string(j);
            std::string refactoredDataFileName = "refactored.tier."+std::to_string(i)+".data."+std::to_string(j)+".bp";
            std::string dataPath = fullDataPath + refactoredDataFileName;
            adios2::Engine data_writer_engine =
                writer_io.Open(dataPath, adios2::Mode::Write); 
            data_writer_engines[idxStr] = data_writer_engine;
        }
        for (size_t j = 0; j < storageTiersECParam_m[i]; j++)
        {
            std::string idxStr = std::to_string(i)+"_"+std::to_string(j);
            std::string refactoredParityFileName = "refactored.tier."+std::to_string(i)+".parity."+std::to_string(j)+".bp";
            std::string parityPath = fullDataPath + refactoredParityFileName;
            adios2::Engine parity_writer_engine =
                writer_io.Open(parityPath, adios2::Mode::Write); 
            parity_writer_engines[idxStr] = parity_writer_engine;
        }        
    }

    std::string fullMetadataPath = storageTiersPaths[0];
    std::string refactoredMetadataFileName = "refactored.md.bp";
    if (!fullMetadataPath.empty() && fullMetadataPath.back() != '/')
    {
        fullMetadataPath += '/';
    }
    fullMetadataPath += refactoredMetadataFileName;
    adios2::Engine metadata_writer_engine =
        writer_io.Open(fullMetadataPath, adios2::Mode::Write);         

    adios2::IO reader_io = adios.DeclareIO("ReaderIO");
    adios2::Engine reader_engine =
           reader_io.Open(inputFileName, adios2::Mode::Read);
    const std::map<std::string, adios2::Params> allVariables =
           reader_io.AvailableVariables();

    for (const auto variablePair : allVariables)
    {
        std::string variableName;
        std::string variableType;
        variableName = variablePair.first;
        variableType = variablePair.second.at("Type");    
        std::vector<uint32_t> storageTiersSizes;
        
        if (variableType == "float")
        {
            auto variable = reader_io.InquireVariable<float>(variableName);
            size_t spaceDimensions = variable.Shape().size();

            std::cout << "read " << spaceDimensions << "D variable " << variableName << std::endl;

            size_t variableSize = 1;
            for (size_t i = 0; i < variable.Shape().size(); i++)
            {
                variableSize *= variable.Shape()[i];
            }

            
            std::vector<float> variableData(variableSize);
            reader_engine.Get(variable, variableData.data(), adios2::Mode::Sync);

            using T = float;
            using T_stream = uint32_t;

            auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
            // auto decomposer = MDR::MGARDHierarchicalDecomposer<T>();
            auto interleaver = MDR::DirectInterleaver<T>();
            // auto interleaver = MDR::SFCInterleaver<T>();
            // auto interleaver = MDR::BlockedInterleaver<T>();
            // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
            auto encoder = MDR::NegaBinaryBPEncoder<T, T_stream>();
            // auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
            // auto compressor = MDR::DefaultLevelCompressor();
            auto compressor = MDR::AdaptiveLevelCompressor(32);
            // auto compressor = MDR::NullLevelCompressor();

            std::vector<uint32_t> dimensions(variable.Shape().size());
            std::vector<T> level_error_bounds;
            std::vector<std::vector<uint8_t*>> level_components;
            std::vector<std::vector<uint32_t>> level_sizes;
            std::vector<std::vector<double>> level_squared_errors;
            std::vector<uint8_t> stopping_indices;

            for (size_t i = 0; i < variable.Shape().size(); i++)
            {
                dimensions[i] = variable.Shape()[i];
                //std::cout << dimensions[i] << " ";
            }
            //std::cout << std::endl;
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
            std::vector<std::vector<double>> level_abs_errors;
            std::vector<std::vector<double>> level_errors;
            MDR::MaxErrorCollector<T> collector = MDR::MaxErrorCollector<T>();
            for(size_t j = 0; j <= target_level; j++)
            {
                auto collected_error = collector.collect_level_error(NULL, 0, level_squared_errors[j].size(), level_error_bounds[j]);
                level_abs_errors.push_back(collected_error);
            }
            level_errors = level_abs_errors;
            auto estimator = MDR::MaxErrorEstimatorOB<T>(spaceDimensions);
            std::vector<uint8_t> idx(target_level+1, 0);
            auto retrieve_order = calculate_retrieve_order(level_sizes, level_errors, 0.0000001, idx, estimator);
            uint32_t totalSize = 0;
            for (size_t i = 0; i < level_sizes.size(); i++)
            {
                for (size_t j = 0; j < level_sizes[i].size(); j++)
                {
                    totalSize += level_sizes[i][j];
                }
            }
            std::cout << totalSize << std::endl;
            uint32_t preSize = 0;
            for (size_t i = 0; i < storageTiersPercentages.size(); i++)
            {
                if (i == storageTiersPercentages.size()-1)
                {
                    storageTiersSizes.push_back(totalSize-preSize);
                }
                else
                {
                    uint32_t oneTierSize = static_cast<uint32_t>(storageTiersPercentages[i]*totalSize);
                    storageTiersSizes.push_back(oneTierSize);
                    preSize += oneTierSize;
                }
            }
            for (size_t i = 0; i < storageTiersSizes.size(); i++)
            {
                std::cout << "tier " << i << ": " << storageTiersSizes[i] << std::endl;
            }
            std::vector<std::vector<uint8_t>> storageTiersValues(storageTiersPaths.size());
            uint32_t currentTier = 0;
            uint32_t currentTierCopidedSize = 0;
            uint32_t totalCopiedSize = 0;
            std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> queryTable;
            for (size_t i = 0; i < retrieve_order.size(); i++)
            {
                uint32_t lid = std::get<0>(retrieve_order[i]);
                uint32_t pid = std::get<1>(retrieve_order[i]);
                std::vector<uint8_t> onePieceValues(level_components[lid][pid], level_components[lid][pid]+level_sizes[lid][pid]);
                // if (i == 0)
                // {
                //     std::cout <<  lid << ", " << pid << ": ";
                //     for (size_t j = 0; j < onePieceValues.size(); j++)
                //     {
                //     std::cout << +onePieceValues[j] << " ";
                //     }
                //     std::cout << std::endl;
                // }
                if (currentTierCopidedSize+level_sizes[lid][pid] <= storageTiersSizes[currentTier])
                {
                    std::cout << "copy " << level_sizes[lid][pid] << " to buffer for tier " << currentTier;
                    storageTiersValues[currentTier].insert(storageTiersValues[currentTier].end(), onePieceValues.begin(), onePieceValues.end());
                    queryTable.push_back(std::make_tuple(lid, pid, currentTier, currentTierCopidedSize, level_sizes[lid][pid]));
                }
                else
                {
                    if (currentTier < storageTiersPaths.size()-1)
                    {
                        std::cout << "tier " << currentTier << "'s buffer is full! start copying to buffer for next tier" << std::endl;
                        currentTier += 1;
                        currentTierCopidedSize = 0;
                    }
                    std::cout << "copy " << level_sizes[lid][pid] << " to buffer for tier " << currentTier;
                    storageTiersValues[currentTier].insert(storageTiersValues[currentTier].end(), onePieceValues.begin(), onePieceValues.end());
                    queryTable.push_back(std::make_tuple(lid, pid, currentTier, currentTierCopidedSize, level_sizes[lid][pid]));
                }
                currentTierCopidedSize += level_sizes[lid][pid];
                totalCopiedSize += level_sizes[lid][pid];
                std::cout << " (" << totalCopiedSize << ")" << std::endl;
            }
            std::vector<uint32_t> queryTableContent;
            for (size_t i = 0; i < queryTable.size(); i++)
            {
                //std::cout << std::get<0>(queryTable[i]) << ", " << std::get<1>(queryTable[i]) << ", " << std::get<2>(queryTable[i]) << ", " << std::get<3>(queryTable[i]) << ", " << std::get<4>(queryTable[i]) << std::endl;
                queryTableContent.push_back(std::get<0>(queryTable[i]));
                queryTableContent.push_back(std::get<1>(queryTable[i]));
                queryTableContent.push_back(std::get<2>(queryTable[i]));
                queryTableContent.push_back(std::get<3>(queryTable[i]));
                queryTableContent.push_back(std::get<4>(queryTable[i]));
            }
            std::string varQueryTableName = variableName+":QueryTable"; 
            adios2::Variable<uint32_t> varQueryTable = writer_io.DefineVariable<uint32_t>(varQueryTableName, {queryTable.size(), 5}, {0, 0}, {queryTable.size(), 5});
            metadata_writer_engine.Put(varQueryTable, queryTableContent.data(), adios2::Mode::Sync);   

            // for (size_t i = 0; i < level_components.size(); i++)
            // {
            //     uint32_t concated_level_size = 0;
            //     for(size_t j = 0; j < level_components[i].size(); j++)
            //     {
            //         concated_level_size += level_sizes[i][j];
            //     }
            //     uint8_t * concated_level_data = (uint8_t *) malloc(concated_level_size);
            //     uint8_t * concated_level_data_pos = concated_level_data;
            //     for(size_t j = 0; j < level_components[i].size(); j++)
            //     {
            //         // if (i == 0 && j == 0)
            //         // {
            //         //     for (size_t k = 0; k < level_sizes[i][j]; k++)
            //         //     {
            //         //         std::cout << +(*(level_components[i][j]+k)) << " ";
            //         //     }
            //         //     std::cout << std::endl;
            //         // }
                    
            //         memcpy(concated_level_data_pos, level_components[i][j], level_sizes[i][j]);
            //         concated_level_data_pos += level_sizes[i][j];
            //     } 
            //     std::string varLevelValuesName = variableName+":Values:level-"+std::to_string(i); 
            //     int assignedTier = static_cast<int>(static_cast<double>(i)*(static_cast<double>(static_cast<double>(storageTiers)/level_components.size())));
            //     //std::cout << "level " << i << ", " << "tier " << assignedTier << std::endl;
            //     adios2::Variable<uint8_t> varLevelValues = writer_io.DefineVariable<uint8_t>(varLevelValuesName, {concated_level_size}, {0}, {concated_level_size});
            //     data_writer_engines[assignedTier].Put(varLevelValues, concated_level_data, adios2::Mode::Sync);
            //     std::string varLevelSizesName = variableName+":Sizes:level-"+std::to_string(i);
            //     adios2::Variable<uint32_t> varLevelSizes = writer_io.DefineVariable<uint32_t>(varLevelSizesName, {level_sizes[i].size()}, {0}, {level_sizes[i].size()});
            //     metadata_writer_engine.Put(varLevelSizes, level_sizes[i].data(), adios2::Mode::Sync);   
            //     std::string varLevelLocationsName = variableName+":Locations:level-"+std::to_string(i);
            //     adios2::Variable<std::string> varLevelLocations = writer_io.DefineVariable<std::string>(varLevelLocationsName);
            //     std::string fullDataPath = storageTiersPaths[assignedTier];
            //     std::string refactoredDataFileName = "refactored.data.bp";
            //     if (!fullDataPath.empty() && fullDataPath.back() != '/')
            //     {
            //         fullDataPath += '/';
            //     }
            //     fullDataPath += refactoredDataFileName;
            //     metadata_writer_engine.Put(varLevelLocations, fullDataPath, adios2::Mode::Sync);         
            // }
            std::string varDimensionsName = variableName+":Dimensions";
            adios2::Variable<uint32_t> varDimensions = writer_io.DefineVariable<uint32_t>(varDimensionsName, {dimensions.size()}, {0}, {dimensions.size()});
            metadata_writer_engine.Put(varDimensions, dimensions.data(), adios2::Mode::Sync);  
            std::string varTypeName = variableName+":Type";
            adios2::Variable<std::string> varType = writer_io.DefineVariable<std::string>(varTypeName);
            metadata_writer_engine.Put(varType, variableType, adios2::Mode::Sync); 

            std::string varLevelsName = variableName+":Levels";
            adios2::Variable<uint32_t> varLevels = writer_io.DefineVariable<uint32_t>(varLevelsName);
            uint32_t numLevels = level_components.size();
            metadata_writer_engine.Put(varLevels, numLevels, adios2::Mode::Sync);            

            std::string varErrorBoundsName = variableName+":ErrorBounds";
            adios2::Variable<T> varErrorBounds = writer_io.DefineVariable<T>(varErrorBoundsName, {level_error_bounds.size()}, {0}, {level_error_bounds.size()});
            metadata_writer_engine.Put(varErrorBounds, level_error_bounds.data(), adios2::Mode::Sync);

            std::string varStopIndicesName = variableName+":StopIndices";
            adios2::Variable<uint8_t> varStopIndices = writer_io.DefineVariable<uint8_t>(varStopIndicesName, {stopping_indices.size()}, {0}, {stopping_indices.size()});
            metadata_writer_engine.Put(varStopIndices, stopping_indices.data(), adios2::Mode::Sync);

            std::vector<double> all_squared_errors;
            for (size_t i = 0; i < level_squared_errors.size(); i++)
            {
                for (size_t j = 0; j < level_squared_errors[i].size(); j++)
                {
                    all_squared_errors.push_back(level_squared_errors[i][j]);
                }    
            }
            std::string varSquaredErrorsName = variableName+":SquaredErrors";
            adios2::Variable<double> varSquaredErrors = writer_io.DefineVariable<double>(varSquaredErrorsName, {level_squared_errors.size(), level_squared_errors[0].size()}, {0, 0}, {level_squared_errors.size(), level_squared_errors[0].size()});
            metadata_writer_engine.Put(varSquaredErrors, all_squared_errors.data(), adios2::Mode::Sync);

            // for (size_t i = 0; i < level_error_bounds.size(); i++)
            // {
            //     //std::string varLevelSizesName = variableName+":Sizes:level-"+std::to_string(i);
            //     //adios2::Variable<uint32_t> varLevelSizes = writer_io.DefineVariable<uint32_t>(varLevelSizesName, {level_sizes[i].size()}, {0}, {level_sizes[i].size()});
            //     //metadata_writer_engine.Put(varLevelSizes, level_sizes[i].data(), adios2::Mode::Sync);  
            //     std::string varLevelErrorBoundsName = variableName+":ErrorBounds:level-"+std::to_string(i);
            //     adios2::Variable<T> varLevelErrorBounds = writer_io.DefineVariable<T>(varLevelErrorBoundsName);
            //     metadata_writer_engine.Put(varLevelErrorBounds, level_error_bounds[i], adios2::Mode::Sync);
            //     std::string varLevelSquaredErrorsName = variableName+":SquaredErrors:level-"+std::to_string(i);
            //     adios2::Variable<double> varLevelSquaredErrors = writer_io.DefineVariable<double>(varLevelSquaredErrorsName, {level_squared_errors[i].size()}, {0}, {level_squared_errors[i].size()});
            //     metadata_writer_engine.Put(varLevelSquaredErrors, level_squared_errors[i].data(), adios2::Mode::Sync);
            // }

            std::string varTiersName = variableName+":Tiers";
            adios2::Variable<uint32_t> varTiers = writer_io.DefineVariable<uint32_t>(varTiersName);
            uint32_t numTiers = storageTiersValues.size();
            metadata_writer_engine.Put(varTiersName, numTiers, adios2::Mode::Sync);  

            // for (size_t i = 0; i < storageTiersValues.size(); i++)
            // {
            //     std::string varStorageTierValuesName = variableName+":Values:"+std::to_string(i); 
            //     adios2::Variable<uint8_t> varStorageTierValues = writer_io.DefineVariable<uint8_t>(varStorageTierValuesName, {storageTiersValues[i].size()}, {0}, {storageTiersValues[i].size()});
            //     data_writer_engines[i].Put(varStorageTierValues, storageTiersValues[i].data(), adios2::Mode::Sync);
            //     std::string varStorageTierLocationsName = variableName+":Tier:"+std::to_string(i);
            //     adios2::Variable<std::string> varStorageTierLocations = writer_io.DefineVariable<std::string>(varStorageTierLocationsName);
            //     std::string fullDataPath = storageTiersPaths[i];
            //     std::string refactoredDataFileName = "refactored.data.bp";
            //     if (!fullDataPath.empty() && fullDataPath.back() != '/')
            //     {
            //         fullDataPath += '/';
            //     }
            //     fullDataPath += refactoredDataFileName;
            //     metadata_writer_engine.Put(varStorageTierLocations, fullDataPath, adios2::Mode::Sync); 
            // }
            for (size_t i = 0; i < storageTiersValues.size(); i++)
            {
                struct ec_args args = {
                    .k = storageTiersECParam_k[i],
                    .m = storageTiersECParam_m[i],
                    .w = storageTiersECParam_w[i],
                    .hd = storageTiersECParam_m[i]+1,
                    .ct = CHKSUM_NONE,
                };
                std::string varECParam_k_Name = variableName+":Tier:"+std::to_string(i)+":K";
                adios2::Variable<int> varECParam_k = writer_io.DefineVariable<int>(varECParam_k_Name); 
                std::string varECParam_m_Name = variableName+":Tier:"+std::to_string(i)+":M";
                adios2::Variable<int> varECParam_m = writer_io.DefineVariable<int>(varECParam_m_Name);
                std::string varECParam_w_Name = variableName+":Tier:"+std::to_string(i)+":W";
                adios2::Variable<int> varECParam_w = writer_io.DefineVariable<int>(varECParam_w_Name);
                std::string varECParam_hd_Name = variableName+":Tier:"+std::to_string(i)+":HD";
                adios2::Variable<int> varECParam_hd = writer_io.DefineVariable<int>(varECParam_hd_Name);

                metadata_writer_engine.Put(varECParam_k, storageTiersECParam_k[i], adios2::Mode::Sync);
                metadata_writer_engine.Put(varECParam_m, storageTiersECParam_m[i], adios2::Mode::Sync);
                metadata_writer_engine.Put(varECParam_w, storageTiersECParam_w[i], adios2::Mode::Sync);
                metadata_writer_engine.Put(varECParam_hd, storageTiersECParam_m[i]+1, adios2::Mode::Sync);

                std::string varECBackendName = variableName+":Tier:"+std::to_string(i)+":ECBackendName";
                adios2::Variable<std::string> varECBackend = writer_io.DefineVariable<std::string>(varECBackendName);
                metadata_writer_engine.Put(varECBackend, ECBackendName, adios2::Mode::Sync);

                int desc = -1;
                int rc = 0;
                //std::cout << "backendID: " << backendID << std::endl;
                desc = liberasurecode_instance_create(backendID, &args);
                //std::cout << desc << std::endl;
                if (-EBACKENDNOTAVAIL == desc) 
                {
                    std::cerr << "backend library not available!" << std::endl;
                    return 1;
                } else if ((args.k + args.m) > EC_MAX_FRAGMENTS) 
                {
                    assert(-EINVALIDPARAMS == desc);
                    std::cerr << "invalid parameters!" << std::endl;
                    return 1;
                } else
                {
                    assert(desc > 0);
                }   
                char **encoded_data = NULL, **encoded_parity = NULL;
                uint64_t encoded_fragment_len = 0;
                char *orig_data = NULL;
                int orig_data_size = storageTiersValues[i].size();
                orig_data = static_cast<char*>(static_cast<void *>(storageTiersValues[i].data()));
                // for (size_t j = 0; j < 10; j++)
                // {
                //     std::cout << storageTiersValues[i][j] << " ";
                // }
                // std::cout << std::endl;
                // for (size_t j = 0; j < 10; j++)
                // {
                //     std::cout << orig_data[j] << " ";
                // }
                // std::cout << std::endl;                
                rc = liberasurecode_encode(desc, orig_data, orig_data_size,
                        &encoded_data, &encoded_parity, &encoded_fragment_len);    
                assert(0 == rc);
                std::cout << "encoded_fragment_len: " << encoded_fragment_len << std::endl;
                std::string varECParam_EncodedFragLen_Name = variableName+":Tier:"+std::to_string(i)+":EncodedFragmentLength";
                adios2::Variable<uint64_t> varECParam_EncodedFragLen = writer_io.DefineVariable<uint64_t>(varECParam_EncodedFragLen_Name);
                metadata_writer_engine.Put(varECParam_EncodedFragLen, encoded_fragment_len, adios2::Mode::Sync);

                size_t frag_header_size =  sizeof(fragment_header_t);
                for (size_t j = 0; j < storageTiersECParam_k[i]; j++)
                {
                    char *frag = NULL;
                    frag = encoded_data[j];
                    assert(frag != NULL);
                    fragment_header_t *header = (fragment_header_t*)frag;
                    assert(header != NULL);

                    fragment_metadata_t metadata = header->meta;
                    assert(metadata.idx == j);
                    assert(metadata.size == encoded_fragment_len - frag_header_size - metadata.frag_backend_metadata_size);
                    assert(metadata.orig_data_size == orig_data_size);
                    assert(metadata.backend_id == backendID);
                    assert(metadata.chksum_mismatch == 0);     

                    std::string varDataValuesName = variableName+":Tier:"+std::to_string(i)+":Data:"+std::to_string(j);   
                    adios2::Variable<char> varDataValues = writer_io.DefineVariable<char>(varDataValuesName, {encoded_fragment_len}, {0}, {encoded_fragment_len});
                    std::string idxStr = std::to_string(i)+"_"+std::to_string(j);
                    data_writer_engines[idxStr].Put(varDataValues, frag, adios2::Mode::Sync);
                    std::string varDataLocationName = variableName+":Tier:"+std::to_string(i)+":Data:"+std::to_string(j)+":Location";
                    adios2::Variable<std::string> varDataLocation = writer_io.DefineVariable<std::string>(varDataLocationName); 
                    std::string fullDataPath = storageTiersPaths[i];
                    if (!fullDataPath.empty() && fullDataPath.back() != '/')
                    {
                        fullDataPath += '/';
                    }
                    std::string refactoredDataFileName = "refactored.tier."+std::to_string(i)+".data."+std::to_string(j)+".bp";
                    std::string dataPath = fullDataPath + refactoredDataFileName;    
                    metadata_writer_engine.Put(varDataLocation, dataPath, adios2::Mode::Sync);                 
                }
                for (size_t j = 0; j < storageTiersECParam_m[i]; j++)
                {
                    char *frag = NULL;
                    frag = encoded_parity[j];
                    assert(frag != NULL);
                    fragment_header_t *header = (fragment_header_t*)frag;
                    assert(header != NULL);

                    fragment_metadata_t metadata = header->meta;
                    assert(metadata.idx == args.k+j);
                    assert(metadata.size == encoded_fragment_len - frag_header_size - metadata.frag_backend_metadata_size);
                    assert(metadata.orig_data_size == orig_data_size);
                    assert(metadata.backend_id == backendID);
                    assert(metadata.chksum_mismatch == 0);     

                    std::string varParityValuesName = variableName+":Tier:"+std::to_string(i)+":Parity:"+std::to_string(j);        
                    adios2::Variable<char> varParityValues = writer_io.DefineVariable<char>(varParityValuesName, {encoded_fragment_len}, {0}, {encoded_fragment_len});  
                    std::string idxStr = std::to_string(i)+"_"+std::to_string(j);
                    parity_writer_engines[idxStr].Put(varParityValues, frag, adios2::Mode::Sync);     
                    std::string varParityLocationName = variableName+":Tier:"+std::to_string(i)+":Parity:"+std::to_string(j)+":Location";
                    adios2::Variable<std::string> varParityLocation = writer_io.DefineVariable<std::string>(varParityLocationName); 
                    std::string fullParityPath = storageTiersPaths[i];
                    if (!fullParityPath.empty() && fullParityPath.back() != '/')
                    {
                        fullParityPath += '/';
                    }
                    std::string refactoredParityFileName = "refactored.tier."+std::to_string(i)+".parity."+std::to_string(j)+".bp";
                    std::string parityPath = fullParityPath + refactoredParityFileName;    
                    metadata_writer_engine.Put(varParityLocation, parityPath, adios2::Mode::Sync);  

                }
                

                rc = liberasurecode_encode_cleanup(desc, encoded_data, encoded_parity);
                assert(rc == 0);    
                assert(0 == liberasurecode_instance_destroy(desc));                 
            }
            
        } 
    }    
    for (auto it : data_writer_engines)
    {
        it.second.Close();
    }
    for (auto it : parity_writer_engines)
    {
        it.second.Close();
    }
    metadata_writer_engine.Close();
    
    reader_engine.Close();
}