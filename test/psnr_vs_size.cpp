#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <queue>
#include <numeric>
#include <type_traits>

#include <adios2.h>

#include "utils.hpp"
#include "../include/Decomposer/Decomposer.hpp"
#include "../include/Interleaver/Interleaver.hpp"
#include "../include/BitplaneEncoder/BitplaneEncoder.hpp"
#include "../include/Retriever/Retriever.hpp"
#include "../include/ErrorEstimator/ErrorEstimator.hpp"
#include "../include/ErrorCollector/ErrorCollector.hpp"
#include "../include/SizeInterpreter/SizeInterpreter.hpp"
#include "../include/LosslessCompressor/LevelCompressor.hpp"
#include "../include/RefactorUtils.hpp"


#include <erasurecode.h>
#include <erasurecode_helpers.h>
#include <config_liberasurecode.h>
#include <erasurecode_stdinc.h>
#include <erasurecode_version.h>

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

std::vector<std::vector<uint32_t>> get_level_sizes(uint32_t levels, const std::vector<std::vector<uint32_t>>& query_table)
{
    std::vector<std::vector<uint32_t>> level_sizes(levels);
    for (size_t i = 0; i < query_table.size(); i++)
    {
        level_sizes[query_table[i][0]].push_back(query_table[i][4]);
    }
    return level_sizes;   
}

std::vector<uint32_t> get_tiers_retrieve_info(uint32_t tiers, const std::vector<std::vector<uint32_t>>& query_table, uint32_t &stop_level, uint32_t &stop_plane, const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, double tolerance, std::vector<uint8_t>& index, MDR::MaxErrorEstimatorOB<float> error_estimator) 
{
    std::vector<uint32_t> tiers_retrieve_info(tiers, 0);
    size_t num_levels = level_sizes.size();
    double accumulated_error = 0;
    for(size_t i = 0; i < num_levels; i++){
        accumulated_error += error_estimator.estimate_error(level_errors[i][index[i]], i);
    }
    std::priority_queue<UnitErrorGain, std::vector<UnitErrorGain>, CompareUniteErrorGain> heap;
    // identify minimal level
    double min_error = accumulated_error;
    for(size_t i = 0; i < num_levels; i++){
        min_error -= error_estimator.estimate_error(level_errors[i][index[i]], i);
        min_error += error_estimator.estimate_error(level_errors[i].back(), i);
        // fetch the first component
        if(index[i] == 0){
            //retrieve_sizes[i] += level_sizes[i][index[i]];
            accumulated_error -= error_estimator.estimate_error(level_errors[i][index[i]], i);
            accumulated_error += error_estimator.estimate_error(level_errors[i][index[i] + 1], i);
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
        if(min_error < tolerance){
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
        stop_level = i;
        stop_plane = index[i];
        index[i]++;
        if(index[i] != level_sizes[i].size())
        {
            double error_gain = error_estimator.estimate_error_gain(accumulated_error, level_errors[i][index[i]], level_errors[i][index[i]+1], i);
            heap.push(UnitErrorGain(error_gain/level_sizes[i][index[i]], i));
        }
        //std::cout << i;
        //std::cout << i << ", " << retrieve_sizes[i] << std::endl;
    }
    //std::cout << stop_level << ", " << stop_plane << std::endl;
    for (size_t i = 0; i < query_table.size(); i++)
    {
        tiers_retrieve_info[query_table[i][2]] += query_table[i][4];
        if (query_table[i][0] == stop_level && query_table[i][1] == stop_plane)
        {
            break;
        }
    }
    return tiers_retrieve_info;
}

void shuffle(std::vector<size_t> &arr, size_t n, unsigned int seed)
{
    if (n > 1) 
    {
        size_t i;
        srand(seed);
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          size_t t = arr[j];
          arr[j] = arr[i];
          arr[i] = t;
        }
    }
}

std::vector<size_t> randomly_mark_dev_as_unavailable(size_t total_dev, size_t unavailable_dev, unsigned int seed)
{
    std::vector<size_t> dev_id_list(total_dev);
    std::iota(dev_id_list.begin(), dev_id_list.end(), 0);
    // for (size_t i = 0; i < dev_id_list.size(); i++)
    // {
    //     std::cout << dev_id_list[i] << " ";
    // }
    // std::cout << std::endl;
    shuffle(dev_id_list, total_dev, seed);
    // for (size_t i = 0; i < dev_id_list.size(); i++)
    // {
    //     std::cout << dev_id_list[i] << " ";
    // }
    // std::cout << std::endl;
    std::vector<size_t> unavailable_dev_list(unavailable_dev);
    for (size_t i = 0; i < unavailable_dev; i++)
    {
        unavailable_dev_list[i] = dev_id_list[i];
    }
    
    return unavailable_dev_list;
}

int main(int argc, char *argv[])
{
    std::string inputFileName;
    std::string variableName;
    int error_mode = 0;
    //int num_tolerance = 0;
    //std::vector<double> tolerance;
    // int totalDevices = 0;
    // int unavaialbleDevices = 0; 
    double s;
    std::string rawDataFileName;
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
        else if (arg == "-var" || arg == "--variable")
        {
            if (i+1 < argc)
            {
                variableName = argv[i+1];
            }
            else
            {
                std::cerr << "--variable option requires one argument." << std::endl;
                return 1;
            } 
        }
        else if (arg == "-em" || arg == "--errormode")
        {
            if (i+1 < argc)
            {
                error_mode = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--errormode option requires one argument." << std::endl;
                return 1;
            }            
        }  
        // else if (arg == "-t" || arg == "--totaldev")
        // {
        //     if (i+1 < argc)
        //     {
        //         totalDevices = atoi(argv[i+1]);
        //         if (totalDevices < 0)
        //         {
        //             std::cerr << "--totaldev option must be greater than 0." << std::endl;
        //             return 1;
        //         }
                
        //     }
        //     else
        //     {
        //         std::cerr << "--totaldev option requires one argument." << std::endl;
        //         return 1;
        //     }            
        // } 
        // else if (arg == "-u" || arg == "--unavaldev")
        // {
        //     if (i+1 < argc)
        //     {
        //         unavaialbleDevices = atoi(argv[i+1]);
        //     }
        //     else
        //     {
        //         std::cerr << "--unavaldev option requires one argument." << std::endl;
        //         return 1;
        //     }            
        // } 
        // else if (arg == "-tlrnc" || arg == "--tolerance")
        // {
        //     if (i+1 < argc)
        //     {
        //         num_tolerance = atoi(argv[i+1]);
        //     }
        //     else
        //     {
        //         std::cerr << "--tolerance option requires [# of tolerance] to be set first." << std::endl;
        //         return 1;
        //     }            
        //     if (num_tolerance)
        //     {
        //         if (i+1+num_tolerance < argc)
        //         {
        //             for (size_t j = i+2; j < i+2+num_tolerance; j++)
        //             {
        //                 tolerance.push_back(atof(argv[j]));
        //             }  
        //         }
        //         else
        //         {
        //             std::cerr << "--tolerance option requires [# of tolerance] arguments." << std::endl;
        //             return 1;
        //         } 
        //     }
        //     else
        //     {
        //         std::cerr << "--[# of tolerance] needs to be greater than zero." << std::endl;
        //         return 1;
        //     }
        // }
        else if (arg == "-s")
        {
            if (i+1 < argc)
            {
                s = atof(argv[i+1]);
            }
            else
            {
                std::cerr << "-s option requires one argument." << std::endl;
                return 1;
            }            
        }    
        else if (arg == "-r" || arg == "--rawdata")
        {
            if (i+1 < argc)
            {
                rawDataFileName = argv[i+1];
            }
            else
            {
                std::cerr << "--rawdata option requires one argument." << std::endl;
                return 1;
            }            
        }     
    }

    adios2::ADIOS adios;
    adios2::IO reader_io = adios.DeclareIO("ReaderIO");
    adios2::Engine metadata_reader_engine =
           reader_io.Open(inputFileName, adios2::Mode::Read);  
    std::cout << variableName+":Dimensions" << std::endl;
    auto varDimensions = reader_io.InquireVariable<uint32_t>(variableName+":Dimensions");
    std::vector<uint32_t> dimensions(varDimensions.Shape()[0]);
    metadata_reader_engine.Get(varDimensions, dimensions.data(), adios2::Mode::Sync);
    std::cout << "dimensions: ";
    for (size_t i = 0; i < dimensions.size(); i++)
    {
        std::cout << dimensions[i] << " ";
    }
    std::cout << std::endl;
    uint32_t levels;
    auto varLevels = reader_io.InquireVariable<uint32_t>(variableName+":Levels");
    metadata_reader_engine.Get(varLevels, levels, adios2::Mode::Sync);
    std::cout << "# of levels: " << levels << std::endl;

    uint32_t tiers;
    auto varTiers = reader_io.InquireVariable<uint32_t>(variableName+":Tiers");
    metadata_reader_engine.Get(varTiers, tiers, adios2::Mode::Sync);
    std::cout << "# of tiers: " << tiers << std::endl;
    // std::vector<std::string> tier_locations(tiers);
    // for (size_t i = 0; i < tiers; i++)
    // {
    //     std::string varOneTierLocationName = variableName+":Tier:"+std::to_string(i);
    //     auto varOneTierLocation = reader_io.InquireVariable<std::string>(varOneTierLocationName);
    //     std::string one_tier_location;
    //     metadata_reader_engine.Get(varOneTierLocation, one_tier_location, adios2::Mode::Sync);
    //     tier_locations[i] = one_tier_location;         
    // }
    std::vector<int32_t> storageTiersECParam_k(tiers);
    std::vector<int32_t> storageTiersECParam_m(tiers);
    std::vector<int32_t> storageTiersECParam_w(tiers);
    std::vector<int32_t> storageTiersECParam_hd(tiers);
    std::vector<std::vector<std::string>> storageTiersDataLocations(tiers);
    std::vector<std::vector<std::string>> storageTiersParityLocations(tiers);
    for (size_t i = 0; i < tiers; i++)
    {
        int32_t oneTierECParam_k;
        std::string varOneTierECParam_k_Name = variableName+":Tier:"+std::to_string(i)+":K";
        auto varOneTierECParam_k = reader_io.InquireVariable<int32_t>(varOneTierECParam_k_Name);
        metadata_reader_engine.Get(varOneTierECParam_k, oneTierECParam_k, adios2::Mode::Sync);
        storageTiersECParam_k[i] = oneTierECParam_k;
        int32_t oneTierECParam_m;
        std::string varOneTierECParam_m_Name = variableName+":Tier:"+std::to_string(i)+":M";
        auto varOneTierECParam_m = reader_io.InquireVariable<int32_t>(varOneTierECParam_m_Name);
        metadata_reader_engine.Get(varOneTierECParam_m, oneTierECParam_m, adios2::Mode::Sync);
        storageTiersECParam_m[i] = oneTierECParam_m;
        int32_t oneTierECParam_w;
        std::string varOneTierECParam_w_Name = variableName+":Tier:"+std::to_string(i)+":W";
        auto varOneTierECParam_w = reader_io.InquireVariable<int32_t>(varOneTierECParam_w_Name);
        metadata_reader_engine.Get(varOneTierECParam_w, oneTierECParam_w, adios2::Mode::Sync);
        storageTiersECParam_w[i] = oneTierECParam_w;
        int32_t oneTierECParam_hd;
        std::string varOneTierECParam_hd_Name = variableName+":Tier:"+std::to_string(i)+":HD";
        auto varOneTierECParam_hd = reader_io.InquireVariable<int32_t>(varOneTierECParam_hd_Name);
        metadata_reader_engine.Get(varOneTierECParam_hd, oneTierECParam_hd, adios2::Mode::Sync);
        storageTiersECParam_hd[i] = oneTierECParam_hd;
        std::cout << "tier " << i << ": " << storageTiersECParam_k[i] << "(k), " << storageTiersECParam_m[i] << "(m), " << storageTiersECParam_w[i] << "(w), " << storageTiersECParam_hd[i] << "(hd)" << std::endl;

        for (size_t j = 0; j < storageTiersECParam_k[i]; j++)
        {
            std::string varDataLocationName = variableName+":Tier:"+std::to_string(i)+":Data:"+std::to_string(j)+":Location";
            auto varDataLocation = reader_io.InquireVariable<std::string>(varDataLocationName);
            std::string data_location;
            metadata_reader_engine.Get(varDataLocation, data_location, adios2::Mode::Sync);
            storageTiersDataLocations[i].push_back(data_location);
            std::cout << "  location of data " << j << ": " << storageTiersDataLocations[i][j] << std::endl;
        }
        for (size_t j = 0; j < storageTiersECParam_m[i]; j++)
        {
            std::string varParityLocationName = variableName+":Tier:"+std::to_string(i)+":Parity:"+std::to_string(j)+":Location";
            auto varParityLocation = reader_io.InquireVariable<std::string>(varParityLocationName);
            std::string parity_location;
            metadata_reader_engine.Get(varParityLocation, parity_location, adios2::Mode::Sync);
            storageTiersParityLocations[i].push_back(parity_location);
            std::cout << "  location of parity " << j << ": " << storageTiersParityLocations[i][j] << std::endl;
        }
        
    }




    std::string variableType;
    auto varType = reader_io.InquireVariable<std::string>(variableName+":Type");
    metadata_reader_engine.Get(varType, variableType, adios2::Mode::Sync);
    std::cout << "variable type: " << variableType << std::endl;

    auto varQueryTable = reader_io.InquireVariable<uint32_t>(variableName+":QueryTable");
    size_t queryTableSize = 1;
    //std::cout << "query table shape: ";
    for (size_t i = 0; i < varQueryTable.Shape().size(); i++)
    {
        //std::cout << varQueryTable.Shape()[i] << " ";
        queryTableSize *= varQueryTable.Shape()[i];
    }
    //std::cout << std::endl;
    std::vector<uint32_t> queryTableContent(queryTableSize);
    metadata_reader_engine.Get(varQueryTable, queryTableContent.data(), adios2::Mode::Sync);
    std::vector<std::vector<uint32_t>> queryTable(varQueryTable.Shape()[0]);
    for (size_t i = 0; i < varQueryTable.Shape()[0]; i++)
    {
        queryTable[i].insert(queryTable[i].end(), queryTableContent.begin()+i*varQueryTable.Shape()[1], queryTableContent.begin()+i*varQueryTable.Shape()[1]+varQueryTable.Shape()[1]);
        // for (size_t j = 0; j < queryTable[i].size(); j++)
        // {
        //     std::cout << queryTable[i][j] << " ";
        // }
        // std::cout << std::endl;
    }
    std::vector<std::vector<uint32_t>> level_sizes = get_level_sizes(levels, queryTable);
    // for (size_t i = 0; i < level_sizes.size(); i++)
    // {
    //     std::cout << "level " << i << ": ";
    //     for (size_t j = 0; j < level_sizes[i].size(); j++)
    //     {
    //         std::cout << level_sizes[i][j] << " ";
    //     }
    //     std::cout << std::endl;   
    // }
    
    auto varSquaredErrors = reader_io.InquireVariable<double>(variableName+":SquaredErrors");
    size_t squaredErrorsSize = 1;
    //std::cout << "squared errors shape: ";
    for (size_t i = 0; i < varSquaredErrors.Shape().size(); i++)
    {
        //std::cout << varSquaredErrors.Shape()[i] << " ";
        squaredErrorsSize *= varSquaredErrors.Shape()[i];
    }
    //std::cout << std::endl;
    std::vector<double> squaredErrors(squaredErrorsSize);
    metadata_reader_engine.Get(varSquaredErrors, squaredErrors.data(), adios2::Mode::Sync);
    std::vector<std::vector<double>> level_squared_errors(levels);
    if (levels != varSquaredErrors.Shape()[0])
    {
        std::cerr << "# of levels is not consistent in the metadata!" << std::endl;
        return 1;
    }
    size_t pos = 0;
    for (size_t i = 0; i < levels; i++)
    {
        level_squared_errors[i].insert(level_squared_errors[i].end(), squaredErrors.begin()+pos, squaredErrors.begin()+pos+varSquaredErrors.Shape()[1]);
        pos += varSquaredErrors.Shape()[1];
    }
    // for (size_t i = 0; i < level_squared_errors.size(); i++)
    // {
    //     std::cout << "level " << i << " squared errors: ";
    //     for (size_t j = 0; j < level_squared_errors[i].size(); j++)
    //     {
    //         std::cout << level_squared_errors[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    
    std::vector<uint8_t> level_num_bitplanes(levels, 0);



    if (variableType == "float")
    {
        using T = float;
        using T_stream = uint32_t;

        std::vector<T> level_error_bounds(levels);
        auto varErrorBounds = reader_io.InquireVariable<T>(variableName+":ErrorBounds");
        metadata_reader_engine.Get(varErrorBounds, level_error_bounds.data(), adios2::Mode::Sync);
        // for (size_t i = 0; i < level_error_bounds.size(); i++)
        // {
        //     std::cout << "level " << i << " error bounds: "<< level_error_bounds[i] << std::endl;
        // }

        adios2::Engine rawdata_reader_engine =
            reader_io.Open(rawDataFileName, adios2::Mode::Read);   
        auto rawVariable = reader_io.InquireVariable<T>(variableName);
        size_t rawVariableSize = 1;
        for (size_t i = 0; i < rawVariable.Shape().size(); i++)
        {
            rawVariableSize *= rawVariable.Shape()[i];
        }
        //std::cout << "size of raw data is " << rawVariableSize << std::endl;
        std::vector<T> rawVariableData(rawVariableSize);
        rawdata_reader_engine.Get(rawVariable, rawVariableData.data(), adios2::Mode::Sync);
        rawdata_reader_engine.Close();

        auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
        auto compressor = MDR::DefaultLevelCompressor();

        

        switch(error_mode)
        {
            case 1:
            {
                auto estimator = MDR::SNormErrorEstimator<T>(dimensions.size(), levels - 1, s);
                auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);  

                break;          
            }  
            default:
            {
                auto estimator = MDR::MaxErrorEstimatorOB<T>(dimensions.size());
                auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorOB<T>>(estimator); 
                std::vector<uint32_t> prev_tiers_retrieve_info(tiers, 0); 
                int prev_stop_level = -1;
                int prev_stop_plane = -1;
                std::vector<std::vector<uint8_t>> storageTiersValues(tiers);

                // std::vector<size_t> unavailableDevList(unavaialbleDevices);
                // unavailableDevList = randomly_mark_dev_as_unavailable(totalDevices, unavaialbleDevices, 0);
                // for (size_t i = 0; i < unavailableDevList.size(); i++)
                // {
                //     std::cout << unavailableDevList[i] << " ";
                // }
                // std::cout << std::endl;

                // size_t storageTiersRecovered = 0;
                for (size_t i = 0; i < storageTiersValues.size(); i++)
                {
                    // if (storageTiersECParam_m[i] < unavaialbleDevices)
                    // {
                    //     std::cout << "tier " << i << ": " << storageTiersECParam_m[i] <<  " parity chunks are not enough to recover from " << unavaialbleDevices << " unavaialble devices!" << std::endl;
                    //     break;
                    // }

                    struct ec_args args = {
                        .k = storageTiersECParam_k[i],
                        .m = storageTiersECParam_m[i],
                        .w = storageTiersECParam_w[i],
                        .hd = storageTiersECParam_m[i]+1,
                        .ct = CHKSUM_NONE,
                    };

                    uint64_t encoded_fragment_len;
                    std::string varECParam_EncodedFragLen_Name = variableName+":Tier:"+std::to_string(i)+":EncodedFragmentLength";
                    adios2::Variable<uint64_t> varECParam_EncodedFragLen = reader_io.InquireVariable<uint64_t>(varECParam_EncodedFragLen_Name);
                    metadata_reader_engine.Get(varECParam_EncodedFragLen, encoded_fragment_len, adios2::Mode::Sync);

                    ec_backend_id_t backendID;
                    std::string ECBackendName;
                    std::string varECBackendName = variableName+":Tier:"+std::to_string(i)+":ECBackendName";
                    adios2::Variable<std::string> varECBackend = reader_io.InquireVariable<std::string>(varECBackendName);
                    metadata_reader_engine.Get(varECBackend, ECBackendName, adios2::Mode::Sync);
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

                    int rc = 0;
                    int desc = -1;
                    uint64_t decoded_data_len = 0;
                    char *decoded_data = NULL;
                    char **avail_frags = NULL;
                    int num_avail_frags = 0;
                    avail_frags = (char **)malloc((storageTiersECParam_k[i] + storageTiersECParam_m[i]) * sizeof(char *));
                    if (avail_frags == NULL)
                    {
                        num_avail_frags = -1;
                        std::cerr << "memory allocation for avail_frags failed!" << std::endl;
                        return 1;
                    }
                    desc = liberasurecode_instance_create(backendID, &args);
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

                    for (size_t j = 0; j < storageTiersECParam_k[i]; j++)
                    {
                        /* check if data chunks are avaialble */
                        // if (std::find(unavailableDevList.begin(), unavailableDevList.end(), j) != unavailableDevList.end())
                        // {
                        //     std::cout << "data chunk " << j << " is unavailable since device " << j << " failed! skip!" << std::endl;
                        //     continue;
                        // }
                        adios2::Engine data_reader_engine =
                            reader_io.Open(storageTiersDataLocations[i][j], adios2::Mode::Read); 
                        std::string varDataValuesName = variableName+":Tier:"+std::to_string(i)+":Data:"+std::to_string(j);
                        auto varDataValues = reader_io.InquireVariable<char>(varDataValuesName);
                        // for (size_t k = 0; k < varDataValues.Shape().size(); k++)
                        // {
                        //     std::cout << varDataValues.Shape()[k] << " ";
                        // }
                        // std::cout << std::endl;
                        //std::vector<char> encodedValues(varDataValues.Shape()[0]);
                        avail_frags[num_avail_frags] = (char *)malloc(varDataValues.Shape()[0]*sizeof(char));
                        data_reader_engine.Get(varDataValues, avail_frags[num_avail_frags], adios2::Mode::Sync);
                        data_reader_engine.Close();
                        //avail_frags[j] = encodedValues.data();
                        num_avail_frags++;
                    }
                    for (size_t j = 0; j < storageTiersECParam_m[i]; j++)
                    {
                        /* check if parity chunks are avaialble */
                        // if (std::find(unavailableDevList.begin(), unavailableDevList.end(), j+storageTiersECParam_k[i]) != unavailableDevList.end())
                        // {
                        //     std::cout << "parity chunk " << j << " is unavailable since device " << j+storageTiersECParam_k[i] << " failed! skip!" << std::endl;
                        //     continue;
                        // }
                        adios2::Engine parity_reader_engine =
                            reader_io.Open(storageTiersParityLocations[i][j], adios2::Mode::Read); 
                        std::string varParityValuesName = variableName+":Tier:"+std::to_string(i)+":Parity:"+std::to_string(j);
                        auto varParityValues = reader_io.InquireVariable<char>(varParityValuesName);
                        // for (size_t k = 0; k < varParityValues.Shape().size(); k++)
                        // {
                        //     std::cout << varParityValues.Shape()[k] << " ";
                        // }
                        // std::cout << std::endl;
                        //std::vector<char> encodedValues(varParityValues.Shape()[0]);
                        avail_frags[num_avail_frags] = (char *)malloc(varParityValues.Shape()[0]*sizeof(char));
                        parity_reader_engine.Get(varParityValues, avail_frags[num_avail_frags], adios2::Mode::Sync);
                        parity_reader_engine.Close();
                        //avail_frags[j+storageTiersECParam_k[i]] = encodedValues.data();
                        num_avail_frags++;
                    }
                    //std::cout << num_avail_frags << std::endl;
                    assert(num_avail_frags > 0);

                    rc = liberasurecode_decode(desc, avail_frags, num_avail_frags,
                                            encoded_fragment_len, 1,
                                            &decoded_data, &decoded_data_len);   
                    assert(0 == rc);

                    uint8_t *tmp = static_cast<uint8_t*>(static_cast<void *>(decoded_data));  
                    std::vector<uint8_t> oneTierDecodedData(tmp, tmp+decoded_data_len);
                    storageTiersValues[i] = oneTierDecodedData;

                    rc = liberasurecode_decode_cleanup(desc, decoded_data);
                    assert(rc == 0);

                    assert(0 == liberasurecode_instance_destroy(desc));

                    free(avail_frags);

                    // storageTiersRecovered++;
                }
                // std::cout << storageTiersRecovered << " storage tiers recovered!" << std::endl;
                // if (storageTiersRecovered == 0)
                // {
                //     std::cerr << "no storage tier is recovered! all data is unavailable!" << std::endl;
                //     return 1;
                // }
                
                // for (size_t i = 0; i < storageTiersRecovered; i++)
                // {
                //     std::cout << "tier " << i << " size: " << storageTiersValues[i].size() << std::endl;
                //     adios2::Engine tmp_data_reader_engine = 
                //         reader_io.Open("/Users/lwk/Research/Projects/ldrd-esamr/tiers-tmp/tier-"+std::to_string(i)+"/refactored.data.bp", adios2::Mode::Read); 
                //     std::string varTmpName = "variable:Values:"+std::to_string(i);
                //     auto varTmp = reader_io.InquireVariable<uint8_t>(varTmpName);
                //     std::vector<uint8_t> var_tmp_values(varTmp.Shape()[0]);
                //     tmp_data_reader_engine.Get(varTmp, var_tmp_values.data(), adios2::Mode::Sync);
                //     tmp_data_reader_engine.Close();
                //     assert(storageTiersValues[i].size() == varTmp.Shape()[0]);

                //     for (size_t j = 0; j < storageTiersValues[i].size(); j++)
                //     {
                //         // if (j < 10)
                //         // {
                //         //     std::cout << +storageTiersValues[i][j] << ", " << +var_tmp_values[j] << std::endl;
                //         // }
                        
                //         if (storageTiersValues[i][j] != var_tmp_values[j])
                //         {
                //             std::cout << j << "th element is not consistent!" << std::endl;
                //         }
                //     }
                    
                // }
                

                
                
                std::vector<std::vector<const uint8_t*>> level_components(levels);
                size_t retrievedSize = 0;
                for (size_t i = 0; i < queryTable.size(); i++)
                {
                    
                    std::cout << i << ": " << queryTable[i][0] << ", " << queryTable[i][1] << ", " << queryTable[i][2] << ", " << queryTable[i][3] << ", " << queryTable[i][4] << std::endl;
                    
                    uint8_t * buffer = (uint8_t *) malloc(queryTable[i][4]);
                    std::copy(storageTiersValues[queryTable[i][2]].begin()+queryTable[i][3], storageTiersValues[queryTable[i][2]].begin()+queryTable[i][3]+queryTable[i][4], buffer);
                    level_components[queryTable[i][0]].push_back(buffer);
                    level_num_bitplanes[queryTable[i][0]]++;

                    retrievedSize += queryTable[i][4];

                    uint8_t target_level = level_error_bounds.size()-1;
                    std::cout << "target level: " << +target_level << std::endl;
                    int skipped_level = 0;
                    for(size_t j = 0; j <= target_level; j++)
                    {
                        if(level_num_bitplanes[target_level-j] != 0)
                        {
                            skipped_level = j;
                            break;
                        }
                    }
                    std::cout << "skipped level: " << skipped_level << std::endl;
                    target_level -= skipped_level;
                    auto level_dims = MDR::compute_level_dims(dimensions, target_level);
                    auto reconstruct_dimensions = level_dims[target_level];
                    uint32_t num_elements = 1;
                    for(const auto& dim:reconstruct_dimensions)
                    {
                        num_elements *= dim;
                    }
                    //std::cout << num_elements << std::endl;
                    std::vector<T> reconstructedData(num_elements, 0);
                    auto level_elements = MDR::compute_level_elements(level_dims, target_level);
                    std::vector<uint32_t> dims_dummy(reconstruct_dimensions.size(), 0);

                    std::cout << "retrieved size: " << retrievedSize << ", target level: " << +target_level << std::endl; 
                    for(size_t j = 0; j <= target_level; j++)
                    {                      
                        std::cout << j << ", " << level_components[j].size() << ", " << level_sizes[j].size() << ", " << +level_num_bitplanes[j] << std::endl; 
                        compressor.decompress_level(level_components[j], level_sizes[j], 0, level_num_bitplanes[j]);
                        std::cout << "decompress_level" << std::endl;
                        int level_exp = 0;
                        frexp(level_error_bounds[j], &level_exp);
                        std::cout << "frexp" << std::endl;
                        auto level_decoded_data = encoder.progressive_decode(level_components[j], level_elements[j], level_exp, 0, level_num_bitplanes[j], j);
                        std::cout << "progressive_decode" << std::endl;
                        compressor.decompress_release();
                        std::cout << "decompress_release" << std::endl;
                        const std::vector<uint32_t>& prev_dims = (j == 0) ? dims_dummy : level_dims[j-1];
                        interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[j], prev_dims, reconstructedData.data());
                        std::cout << "reposition" << std::endl;
                        free(level_decoded_data);
                        std::cout << "free" << std::endl;
                    }

                    decomposer.recompose(reconstructedData.data(), reconstruct_dimensions, target_level);
                    MGARD::print_statistics(rawVariableData.data(), reconstructedData.data(), rawVariableData.size()); 
                }

                
            }
        }

    }
    metadata_reader_engine.Close();
    
}