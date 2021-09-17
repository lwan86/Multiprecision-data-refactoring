#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include <queue>
#include <unordered_map>

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

std::vector<std::vector<std::tuple<uint32_t, uint32_t>>> get_tiers_retrieve_info(uint32_t tiers, std::unordered_map<std::string, std::vector<uint32_t>>& query_table_dict, std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> &access_order, const std::vector<std::vector<uint32_t>>& level_sizes, const std::vector<std::vector<double>>& level_errors, double tolerance, std::vector<uint8_t>& index, MDR::MaxErrorEstimatorOB<float> error_estimator) 
{
    std::vector<std::vector<std::tuple<uint32_t, uint32_t>>> tiers_retrieve_info(tiers);
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
            std::string dictKey = std::to_string(i)+"+"+std::to_string(index[i]);
            size_t tierID = query_table_dict[dictKey][0];
            if (tiers_retrieve_info[tierID].empty())
            {
                tiers_retrieve_info[tierID].push_back(std::make_tuple(query_table_dict[dictKey][1], query_table_dict[dictKey][2]));
            }
            else
            {
                uint32_t prev_start = std::get<0>(tiers_retrieve_info[tierID].back());
                uint32_t prev_count = std::get<1>(tiers_retrieve_info[tierID].back());
                if (prev_start+prev_count == query_table_dict[dictKey][1])
                {
                    tiers_retrieve_info[tierID].pop_back();
                    tiers_retrieve_info[tierID].push_back(std::make_tuple(prev_start, prev_count+query_table_dict[dictKey][2]));
                }
                else
                {
                    tiers_retrieve_info[tierID].push_back(std::make_tuple(query_table_dict[dictKey][1], query_table_dict[dictKey][2]));
                }
            }
            access_order.push_back(std::make_tuple(i, index[i], query_table_dict[dictKey][0], query_table_dict[dictKey][1], query_table_dict[dictKey][2]));

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
        std::string dictKey = std::to_string(i)+"+"+std::to_string(j);
        size_t tierID = query_table_dict[dictKey][0];
        if (tiers_retrieve_info[tierID].empty())
        {
            tiers_retrieve_info[tierID].push_back(std::make_tuple(query_table_dict[dictKey][1], query_table_dict[dictKey][2]));
        }
        else
        {
            uint32_t prev_start = std::get<0>(tiers_retrieve_info[tierID].back());
            uint32_t prev_count = std::get<1>(tiers_retrieve_info[tierID].back());
            if (prev_start+prev_count == query_table_dict[dictKey][1])
            {
                tiers_retrieve_info[tierID].pop_back();
                tiers_retrieve_info[tierID].push_back(std::make_tuple(prev_start, prev_count+query_table_dict[dictKey][2]));
            }
            else
            {
                tiers_retrieve_info[tierID].push_back(std::make_tuple(query_table_dict[dictKey][1], query_table_dict[dictKey][2]));
            }
        }
        access_order.push_back(std::make_tuple(i, index[i], query_table_dict[dictKey][0], query_table_dict[dictKey][1], query_table_dict[dictKey][2]));

        accumulated_error -= error_estimator.estimate_error(level_errors[i][j], i);
        accumulated_error += error_estimator.estimate_error(level_errors[i][j+1], i);
        if(accumulated_error < tolerance)
        {
            tolerance_met = true;
        }
        //std::cout << i << ", " << +index[i] << std::endl;
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
    // for (size_t i = 0; i < query_table.size(); i++)
    // {
    //     tiers_retrieve_info[query_table[i][2]] += query_table[i][4];
    //     if (query_table[i][0] == stop_level && query_table[i][1] == stop_plane)
    //     {
    //         break;
    //     }
    // }
    return tiers_retrieve_info;
}

int main(int argc, char *argv[])
{
    std::string inputFileName;
    std::string variableName;
    int error_mode = 0;
    int num_tolerance = 0;
    std::vector<double> tolerance;
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
        else if (arg == "-tlrnc" || arg == "--tolerance")
        {
            if (i+1 < argc)
            {
                num_tolerance = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--tolerance option requires [# of tolerance] to be set first." << std::endl;
                return 1;
            }            
            if (num_tolerance)
            {
                if (i+1+num_tolerance < argc)
                {
                    for (size_t j = i+2; j < i+2+num_tolerance; j++)
                    {
                        tolerance.push_back(atof(argv[j]));
                    }  
                }
                else
                {
                    std::cerr << "--tolerance option requires [# of tolerance] arguments." << std::endl;
                    return 1;
                } 
            }
            else
            {
                std::cerr << "--[# of tolerance] needs to be greater than zero." << std::endl;
                return 1;
            }
        }
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
    std::vector<std::string> tier_locations(tiers);
    for (size_t i = 0; i < tiers; i++)
    {
        std::string varOneTierLocationName = variableName+":Tier:"+std::to_string(i);
        auto varOneTierLocation = reader_io.InquireVariable<std::string>(varOneTierLocationName);
        std::string one_tier_location;
        metadata_reader_engine.Get(varOneTierLocation, one_tier_location, adios2::Mode::Sync);
        tier_locations[i] = one_tier_location;         
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
    std::unordered_map<std::string, std::vector<uint32_t>> queryTableDictionary;
    for (size_t i = 0; i < varQueryTable.Shape()[0]; i++)
    {
        queryTable[i].insert(queryTable[i].end(), queryTableContent.begin()+i*varQueryTable.Shape()[1], queryTableContent.begin()+i*varQueryTable.Shape()[1]+varQueryTable.Shape()[1]);
        
        // std::cout << i << ": "; 
        // for (size_t j = 0; j < queryTable[i].size(); j++)
        // {
        //     std::cout << queryTable[i][j] << " ";
        // }
        // std::cout << std::endl;
        std::vector<uint32_t> dictVal;
        dictVal.push_back(queryTable[i][2]);
        dictVal.push_back(queryTable[i][3]);
        dictVal.push_back(queryTable[i][4]);
        std::string dictKey = std::to_string(queryTable[i][0])+"+"+std::to_string(queryTable[i][1]);
        queryTableDictionary[dictKey] = dictVal;
        // std::cout << dictKey << ": ";
        // for (size_t j = 0; j < queryTableDictionary[dictKey].size(); j++)
        // {
        //     std::cout << queryTableDictionary[dictKey][j] << " ";
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

    std::vector<uint8_t> stopping_indices(levels);
    auto varStopIndices = reader_io.InquireVariable<uint8_t>(variableName+":StopIndices");
    metadata_reader_engine.Get(varStopIndices, stopping_indices.data(), adios2::Mode::Sync);    


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

        std::vector<T> reconstructedData;

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

                std::vector<std::vector<std::tuple<uint32_t, uint32_t, std::vector<uint8_t>>>> storageTiersValues(tiers);
                
                for (size_t i = 0; i < tolerance.size(); i++)
                {
                    std::cout << "tolerance target: " << tolerance[i] << std::endl;
                    std::vector<T> exist_data(reconstructedData);
                    std::vector<std::vector<double>> level_abs_errors;
                    uint8_t target_level = level_error_bounds.size()-1;
                    MDR::MaxErrorCollector<T> collector = MDR::MaxErrorCollector<T>();
                    std::vector<std::vector<const uint8_t*>> level_components(levels);
                    for(size_t j = 0; j <= target_level; j++)
                    {
                        auto collected_error = collector.collect_level_error(NULL, 0, level_squared_errors[j].size(), level_error_bounds[j]);
                        level_abs_errors.push_back(collected_error);
                    }
                    std::vector<std::vector<double>> level_errors = level_abs_errors;
                    auto prev_level_num_bitplanes(level_num_bitplanes);  

                    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>> access_order;
                    std::vector<std::vector<std::tuple<uint32_t, uint32_t>>> tiers_retrieve_info = get_tiers_retrieve_info(tiers, queryTableDictionary, access_order, level_sizes, level_errors, tolerance[i], level_num_bitplanes, estimator);
                    for (size_t j = 0; j < tiers_retrieve_info.size(); j++)
                    {
                        //std::cout << "tier " << j << ": " << std::endl;
                        if (tiers_retrieve_info[j].empty())
                        {
                            std::cout << "no need to retrieve from tier " << j << std::endl;
                            continue;
                        }
                        
                        for (size_t k = 0; k < tiers_retrieve_info[j].size(); k++)
                        {
                            adios2::Engine one_tier_values_reader_engine =
                                reader_io.Open(tier_locations[j], adios2::Mode::Read);
                            std::string varOneTierValuesName = variableName+":Values:"+std::to_string(j);
                            auto varOneTierValues = reader_io.InquireVariable<uint8_t>(varOneTierValuesName);
                            varOneTierValues.SetSelection({{std::get<0>(tiers_retrieve_info[j][k])}, {std::get<1>(tiers_retrieve_info[j][k])}});
                            std::vector<uint8_t> one_tier_values(std::get<1>(tiers_retrieve_info[j][k]));
                            one_tier_values_reader_engine.Get(varOneTierValues, one_tier_values.data(), adios2::Mode::Sync);
                            std::cout << "retrieve " << std::get<1>(tiers_retrieve_info[j][k]) << " from tier " << j << " (" << tier_locations[j] << ")" << std::endl;
                            storageTiersValues[j].push_back(std::make_tuple(std::get<0>(tiers_retrieve_info[j][k]), std::get<1>(tiers_retrieve_info[j][k]), one_tier_values));
                            one_tier_values_reader_engine.Close();
                        }
                        
                        // std::cout << "tier " << j << ": ";
                        // for (size_t k = 0; k < tiers_retrieve_info[j].size(); k++)
                        // {
                        //     std::cout << std::get<0>(tiers_retrieve_info[j][k]) << " " << std::get<1>(tiers_retrieve_info[j][k]) << ", ";
                        // }
                        // std::cout << std::endl;
                    }
                    
                    for (size_t j = 0; j < access_order.size(); j++)
                    {
                        //std::cout << j << ": " << std::get<0>(access_order[j]) << ", " << std::get<1>(access_order[j]) << ", " << std::get<2>(access_order[j]) << ", " << std::get<3>(access_order[j]) << ", " << std::get<4>(access_order[j]) << std::endl;
                        uint8_t * buffer = (uint8_t *) malloc(std::get<4>(access_order[j]));
                        for (size_t k = 0; k < storageTiersValues[std::get<2>(access_order[j])].size(); k++)
                        {
                            std::tuple<uint32_t, uint32_t, std::vector<uint8_t>> tier_values = storageTiersValues[std::get<2>(access_order[j])][k];
                            //std::cout << "  " << std::get<0>(tier_values) << ", " << std::get<1>(tier_values) << ", " << std::get<3>(access_order[j]) << ", " << std::get<4>(access_order[j]) << std::endl;
                            if (std::get<0>(tier_values) <= std::get<3>(access_order[j]) && std::get<0>(tier_values)+std::get<1>(tier_values) >= std::get<3>(access_order[j])+std::get<4>(access_order[j]))
                            {
                                std::copy(std::get<2>(tier_values).begin()+std::get<3>(access_order[j])-std::get<0>(tier_values), std::get<2>(tier_values).begin()+std::get<3>(access_order[j])-std::get<0>(tier_values)+std::get<4>(access_order[j]), buffer);
                                level_components[std::get<0>(access_order[j])].push_back(buffer);
                            }
                        }
                    }
                    
                    // check whether to reconstruct to full resolution
                    int skipped_level = 0;
                    for(size_t j = 0; j <= target_level; j++)
                    {
                        if(level_num_bitplanes[target_level-j] != 0)
                        {
                            skipped_level = j;
                            break;
                        }
                    }
                    target_level -= skipped_level;
                    auto level_dims = MDR::compute_level_dims(dimensions, target_level);
                    auto reconstruct_dimensions = level_dims[target_level];
                    uint32_t num_elements = 1;
                    for(const auto& dim:reconstruct_dimensions)
                    {
                        num_elements *= dim;
                    }
                    //std::cout << num_elements << std::endl;
                    reconstructedData.clear();
                    reconstructedData = std::vector<T>(num_elements, 0);
                    
                    auto level_elements = MDR::compute_level_elements(level_dims, target_level);

                    // for (size_t j = 0; j < level_elements.size(); j++)
                    // {
                    //     std::cout << level_elements[j] << " ";
                    // }
                    // std::cout << std::endl;
                    
                    std::vector<uint32_t> dims_dummy(reconstruct_dimensions.size(), 0);
                    //std::cout << level_components.size() << std::endl;
                    
                    for(size_t j = 0; j <= target_level; j++)
                    {
                        //std::cout << "level " << j << " components size: "<< level_components[j].size() << std::endl;
                        //std::cout << "level " << j << " sizes: "<< level_sizes[j].size() << std::endl;
                        //std::cout << "level " << j << " prev_level_num_bitplanes: "<< +prev_level_num_bitplanes[j] << std::endl;
                        //std::cout << "level " << j << " level_num_bitplanes: "<< +level_num_bitplanes[j] << std::endl;
                        // for (size_t k = 0; k < level_components[j].size(); k++)
                        // {
                        //     std::cout << j << ", " << k << ": ";
                        //     for (size_t l = 0; l < 20; l++)
                        //     {
                        //         std::cout << +level_components[j][k][l] << " ";
                        //     }
                        //     std::cout << std::endl;
                        // }
                        //std::cout << j << std::endl;
                        // std::cout << "level " << j << " sizes: ";
                        // for (size_t k = 0; k < level_sizes[j].size(); k++)
                        // {
                        //     std::cout << level_sizes[j][k] << " ";
                        // }
                        // std::cout << std::endl;

                        // std::cout << "level " << j << " components before decompress_level: " << std::endl;
                        // for (size_t k = 0; k < level_components[j].size(); k++)
                        // {
                        //     const uint8_t *ptr = level_components[j][k];
                        //     size_t count = 0;
                        //     while (count < 10)
                        //     {
                        //         std::cout << +(*ptr) << ", ";
                        //         ptr++;
                        //         count++;
                        //     }
                        //     std::cout << std::endl;
                        // }
                        compressor.decompress_level(level_components[j], level_sizes[j], prev_level_num_bitplanes[j], level_num_bitplanes[j] - prev_level_num_bitplanes[j], stopping_indices[j]);
                        int level_exp = 0;
                        frexp(level_error_bounds[j], &level_exp);
                        // std::cout << "level " << j << " components after decompress_level: " << std::endl;
                        // for (size_t k = 0; k < level_components[j].size(); k++)
                        // {
                        //     const uint8_t *ptr = level_components[j][k];
                        //     size_t count = 0;
                        //     while (count < 10)
                        //     {
                        //         std::cout << +(*ptr) << ", ";
                        //         ptr++;
                        //         count++;
                        //     }
                        //     std::cout << std::endl;
                        // }
                        auto level_decoded_data = encoder.progressive_decode(level_components[j], level_elements[j], level_exp, prev_level_num_bitplanes[j], level_num_bitplanes[j] - prev_level_num_bitplanes[j], j);
                        compressor.decompress_release();
                        const std::vector<uint32_t>& prev_dims = (j == 0) ? dims_dummy : level_dims[j-1];
                        interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[j], prev_dims, reconstructedData.data());
                        free(level_decoded_data);
                        //std::cout << " pass" << std::endl;
                    }

                    decomposer.recompose(reconstructedData.data(), reconstruct_dimensions, target_level);
                    if(exist_data.size() == reconstructedData.size())
                    {
                        for(size_t j = 0; j < reconstructedData.size(); j++)
                        {
                            reconstructedData[j] += exist_data[j];
                        }                
                    }
                    else if(exist_data.size())
                    {
                        std::cerr << "Reconstruct size changes, not supported yet." << std::endl;
                        std::cerr << "Sizes before reconstruction: " << exist_data.size() << std::endl;
                        std::cerr << "Sizes after reconstruction: " << reconstructedData.size() << std::endl;
                        return 1;
                    }                    
                    MGARD::print_statistics(rawVariableData.data(), reconstructedData.data(), rawVariableData.size());                     
                    
                    /*
                    uint32_t stop_level, stop_plane;
                    std::vector<uint32_t> tiers_retrieve_info = get_tiers_retrieve_info(tiers, queryTable, stop_level, stop_plane, level_sizes, level_errors, tolerance[i], level_num_bitplanes, estimator);
                    uint32_t current_stop_level = stop_level;
                    uint32_t current_stop_plane = stop_plane;

                    for (size_t j = 0; j < tiers_retrieve_info.size(); j++)
                    {
                        if (tiers_retrieve_info[j]-prev_tiers_retrieve_info[j] == 0)
                        {
                            std::cout << "no need to retrieve from tier " << j << std::endl;
                            continue;
                        }
                        //std::cout << tier_locations[j] << std::endl;
                        adios2::Engine one_tier_values_reader_engine =
                            reader_io.Open(tier_locations[j], adios2::Mode::Read); 
                        std::string varOneTierValuesName = variableName+":Values:"+std::to_string(j);
                        auto varOneTierValues = reader_io.InquireVariable<uint8_t>(varOneTierValuesName);
                        varOneTierValues.SetSelection({{prev_tiers_retrieve_info[j]}, {tiers_retrieve_info[j]-prev_tiers_retrieve_info[j]}});
                        std::vector<uint8_t> one_tier_values(tiers_retrieve_info[j]-prev_tiers_retrieve_info[j]);
                        one_tier_values_reader_engine.Get(varOneTierValues, one_tier_values.data(), adios2::Mode::Sync);
                        std::cout << "retrieve " << tiers_retrieve_info[j]-prev_tiers_retrieve_info[j] << " from tier " << j << " (" << tier_locations[j] << ")" << std::endl;
                        storageTiersValues[j].insert(storageTiersValues[j].end(), one_tier_values.begin(), one_tier_values.end());
                        one_tier_values_reader_engine.Close();

                        // std::cout << "tier " << j << ": ";
                        // for (size_t k = 0; k < 54; k++)
                        // {
                        //     std::cout << +storageTiersValues[j][k] << ", ";
                        // }
                        // std::cout << std::endl;    
                    }
                    std::cout << "prev stop level: " << prev_stop_level << ", prev stop plane: " << prev_stop_plane << std::endl;  
                    std::cout << "current stop level: " << current_stop_level << ", current stop plane:  " << current_stop_plane << std::endl;  
                    size_t queryTableStartPos, queryTableEndPos;
                    
                    for (size_t j = 0; j < queryTable.size(); j++)
                    {
                        //std::cout << prev_stop_level << ", " << prev_stop_plane << std::endl;
                        //std::cout << current_stop_level << ", " << current_stop_plane << std::endl;
                        //std::cout << "j: " << j << ", queryTable[j][0]: "<< queryTable[j][0] << ", queryTable[j][1]: " << queryTable[j][1] << std::endl; 
                        
                        if (prev_stop_level == -1 && prev_stop_plane == -1)
                        {
                            queryTableStartPos = 0;
                        }
                        else if ((queryTable[j][0] == prev_stop_level) && (queryTable[j][1] == prev_stop_plane))
                        {
                            queryTableStartPos = j+1;
                        }
                        if ((queryTable[j][0] == current_stop_level) && (queryTable[j][1] == current_stop_plane))
                        {
                            queryTableEndPos = j;
                            break;
                        }
                    }
                    std::cout << "queryTableStartPos: " << queryTableStartPos << ", queryTableEndPos: " << queryTableEndPos << std::endl;
                    //std::cout <<  queryTable[queryTableStartPos][0] << ", " << queryTable[queryTableStartPos][1] << std::endl;
                    //std::cout <<  queryTable[queryTableEndPos][0] << ", " << queryTable[queryTableEndPos][1] << std::endl;
                    //queryTableStartPos = 0;
                    //queryTableEndPos = 114;
                    for (size_t j = queryTableStartPos; j <= queryTableEndPos; j++)
                    {
                        uint8_t * buffer = (uint8_t *) malloc(queryTable[j][4]);
                        std::copy(storageTiersValues[queryTable[j][2]].begin()+queryTable[j][3], storageTiersValues[queryTable[j][2]].begin()+queryTable[j][3]+queryTable[j][4], buffer);
                        std::cout << j << ": " << queryTable[j][0] << ", " << queryTable[j][1] << ", " << queryTable[j][2] << ", " << queryTable[j][3] << ", " << queryTable[j][5] << std::endl;
                        // for (size_t k = 0; k < 20; k++)
                        // {
                        //     std::cout << +(buffer[k]) << " ";
                        // }
                        // std::cout << std::endl;
                        level_components[queryTable[j][0]].push_back(buffer);
                    }
                    
                    // check whether to reconstruct to full resolution
                    int skipped_level = 0;
                    for(size_t j = 0; j <= target_level; j++)
                    {
                        if(level_num_bitplanes[target_level-j] != 0)
                        {
                            skipped_level = j;
                            break;
                        }
                    }
                    target_level -= skipped_level;
                    auto level_dims = MDR::compute_level_dims(dimensions, target_level);
                    auto reconstruct_dimensions = level_dims[target_level];
                    uint32_t num_elements = 1;
                    for(const auto& dim:reconstruct_dimensions)
                    {
                        num_elements *= dim;
                    }
                    //std::cout << num_elements << std::endl;
                    reconstructedData.clear();
                    reconstructedData = std::vector<T>(num_elements, 0);
                    
                    auto level_elements = MDR::compute_level_elements(level_dims, target_level);

                    // for (size_t j = 0; j < level_elements.size(); j++)
                    // {
                    //     std::cout << level_elements[j] << " ";
                    // }
                    // std::cout << std::endl;
                    
                    std::vector<uint32_t> dims_dummy(reconstruct_dimensions.size(), 0);
                    //std::cout << level_components.size() << std::endl;
                    
                    for(size_t j = 0; j <= target_level; j++)
                    {
                        //std::cout << "level " << j << " components size: "<< level_components[j].size() << std::endl;
                        //std::cout << "level " << j << " sizes: "<< level_sizes[j].size() << std::endl;
                        //std::cout << "level " << j << " prev_level_num_bitplanes: "<< +prev_level_num_bitplanes[j] << std::endl;
                        //std::cout << "level " << j << " level_num_bitplanes: "<< +level_num_bitplanes[j] << std::endl;
                        // for (size_t k = 0; k < level_components[j].size(); k++)
                        // {
                        //     std::cout << j << ", " << k << ": ";
                        //     for (size_t l = 0; l < 20; l++)
                        //     {
                        //         std::cout << +level_components[j][k][l] << " ";
                        //     }
                        //     std::cout << std::endl;
                        // }
                        //std::cout << j << std::endl;
                        std::cout << "level " << j << " sizes: ";
                        for (size_t k = 0; k < level_sizes[j].size(); k++)
                        {
                            std::cout << level_sizes[j][k] << " ";
                        }
                        std::cout << std::endl;

                        std::cout << "level " << j << " components before decompress_level: " << std::endl;
                        for (size_t k = 0; k < level_components[j].size(); k++)
                        {
                            const uint8_t *ptr = level_components[j][k];
                            size_t count = 0;
                            while (count < 10)
                            {
                                std::cout << +(*ptr) << ", ";
                                ptr++;
                                count++;
                            }
                            std::cout << std::endl;
                        }
                        compressor.decompress_level(level_components[j], level_sizes[j], prev_level_num_bitplanes[j], level_num_bitplanes[j] - prev_level_num_bitplanes[j], stopping_indices[j]);
                        int level_exp = 0;
                        frexp(level_error_bounds[j], &level_exp);
                        std::cout << "level " << j << " components after decompress_level: " << std::endl;
                        for (size_t k = 0; k < level_components[j].size(); k++)
                        {
                            const uint8_t *ptr = level_components[j][k];
                            size_t count = 0;
                            while (count < 10)
                            {
                                std::cout << +(*ptr) << ", ";
                                ptr++;
                                count++;
                            }
                            std::cout << std::endl;
                        }
                        auto level_decoded_data = encoder.progressive_decode(level_components[j], level_elements[j], level_exp, prev_level_num_bitplanes[j], level_num_bitplanes[j] - prev_level_num_bitplanes[j], j);
                        compressor.decompress_release();
                        const std::vector<uint32_t>& prev_dims = (j == 0) ? dims_dummy : level_dims[j-1];
                        interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[j], prev_dims, reconstructedData.data());
                        free(level_decoded_data);
                        //std::cout << " pass" << std::endl;
                    }

                    decomposer.recompose(reconstructedData.data(), reconstruct_dimensions, target_level);
                    if(exist_data.size() == reconstructedData.size())
                    {
                        for(size_t j = 0; j < reconstructedData.size(); j++)
                        {
                            reconstructedData[j] += exist_data[j];
                        }                
                    }
                    else if(exist_data.size())
                    {
                        std::cerr << "Reconstruct size changes, not supported yet." << std::endl;
                        std::cerr << "Sizes before reconstruction: " << exist_data.size() << std::endl;
                        std::cerr << "Sizes after reconstruction: " << reconstructedData.size() << std::endl;
                        return 1;
                    }                    
                    MGARD::print_statistics(rawVariableData.data(), reconstructedData.data(), rawVariableData.size()); 
                    prev_tiers_retrieve_info = tiers_retrieve_info;  
                    prev_stop_level = current_stop_level;
                    prev_stop_plane = current_stop_plane;   

                    */             
                }
            }
        }

    }
    metadata_reader_engine.Close();
    
}