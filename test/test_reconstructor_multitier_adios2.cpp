#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>

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
    std::string variableType;
    auto varType = reader_io.InquireVariable<std::string>(variableName+":Type");
    metadata_reader_engine.Get(varType, variableType, adios2::Mode::Sync);
    std::cout << "variable type: " << variableType << std::endl;


    if (variableType == "float")
    {
        using T = float;
        using T_stream = uint32_t;
        std::vector<T> level_error_bounds(levels);
        std::vector<std::vector<double>> level_squared_errors(levels);
        std::vector<std::vector<uint32_t>> level_sizes(levels);
        std::vector<std::string> level_locations(levels);
        std::vector<uint8_t> level_num_bitplanes;
        level_num_bitplanes = std::vector<uint8_t>(levels, 0);
        for (size_t i = 0; i < levels; i++)
        {
            std::string varOneLevelSizesName = variableName+":Sizes:level-"+std::to_string(i);
            auto varOneLevelSizes = reader_io.InquireVariable<uint32_t>(varOneLevelSizesName);
            std::vector<uint32_t> one_level_sizes(varOneLevelSizes.Shape()[0]);
            metadata_reader_engine.Get(varOneLevelSizes, one_level_sizes.data(), adios2::Mode::Sync);
            level_sizes[i] = one_level_sizes;

            std::string varOneLevelSquaredErrorsName = variableName+":SquaredErrors:level-"+std::to_string(i);
            auto varOneLevelSquaredErrors = reader_io.InquireVariable<double>(varOneLevelSquaredErrorsName);
            std::vector<double> one_level_squared_errors(varOneLevelSquaredErrors.Shape()[0]);
            metadata_reader_engine.Get(varOneLevelSquaredErrors, one_level_squared_errors.data(), adios2::Mode::Sync);
            level_squared_errors[i] = one_level_squared_errors;

            std::string varOneLevelErrorBoundsName = variableName+":ErrorBounds:level-"+std::to_string(i);
            auto varOneLevelErrorBounds = reader_io.InquireVariable<T>(varOneLevelErrorBoundsName);
            T one_level_error_bounds;
            metadata_reader_engine.Get(varOneLevelErrorBounds, one_level_error_bounds, adios2::Mode::Sync);
            level_error_bounds[i] = one_level_error_bounds;    

            std::string varOneLevelLocationName = variableName+":Locations:level-"+std::to_string(i);
            auto varOneLevelLocation = reader_io.InquireVariable<std::string>(varOneLevelLocationName);
            std::string one_level_location;
            metadata_reader_engine.Get(varOneLevelLocation, one_level_location, adios2::Mode::Sync);
            level_locations[i] = one_level_location;   
            //std::cout << "level " << i << ": " <<  level_locations[i] << std::endl;             
        }
        adios2::Engine rawdata_reader_engine =
            reader_io.Open(rawDataFileName, adios2::Mode::Read);   
        auto rawVariable = reader_io.InquireVariable<T>(variableName);
        size_t rawVariableSize = 1;
        for (size_t i = 0; i < rawVariable.Shape().size(); i++)
        {
            rawVariableSize *= rawVariable.Shape()[i];
        }
        std::cout << "size of raw data is " << rawVariableSize << std::endl;
        std::vector<T> rawVariableData(rawVariableSize);
        rawdata_reader_engine.Get(rawVariable, rawVariableData.data(), adios2::Mode::Sync);
        rawdata_reader_engine.Close();

        auto decomposer = MDR::MGARDOrthoganalDecomposer<T>();
        auto interleaver = MDR::DirectInterleaver<T>();
        auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
        auto compressor = MDR::DefaultLevelCompressor();

        std::vector<T> reconstructedData;
        struct timespec start, end;
        int err = 0;
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
                std::vector<uint32_t> offsets = std::vector<uint32_t>(levels, 0);
                for (size_t i = 0; i < tolerance.size(); i++)
                {
                    std::cout << "Start reconstruction" << std::endl;
                    err = clock_gettime(CLOCK_REALTIME, &start);

                    std::vector<T> exist_data(reconstructedData);
                    std::vector<std::vector<double>> level_abs_errors;
                    uint8_t target_level = level_error_bounds.size() - 1;
                    std::vector<std::vector<double>>& level_errors = level_squared_errors;
                    MDR::MaxErrorCollector<T> collector = MDR::MaxErrorCollector<T>();
                    for(size_t j = 0; j <= target_level; j++)
                    {
                        auto collected_error = collector.collect_level_error(NULL, 0, level_squared_errors[j].size(), level_error_bounds[j]);
                        level_abs_errors.push_back(collected_error);
                    }
                    level_errors = level_abs_errors;
                    auto prev_level_num_bitplanes(level_num_bitplanes);
                    for (size_t j = 0; j < levels; j++)
                    {
                        for (size_t k = 0; k < level_sizes[j].size(); k++)
                        {
                            std::cout << level_sizes[j][k] << " ";
                        }
                        std::cout << std::endl;
                    }
                    for (size_t j = 0; j < levels; j++)
                    {
                        for (size_t k = 0; k < level_errors[j].size(); k++)
                        {
                            std::cout << level_errors[j][k] << " ";
                        }
                        std::cout << std::endl;
                    }
                    
                    auto retrieve_sizes = interpreter.interpret_retrieve_size(level_sizes, level_errors, tolerance[i], level_num_bitplanes);
                    std::vector<std::vector<const uint8_t*>> level_components;
                    std::vector<uint8_t*> concated_level_components;
                    uint32_t total_retrieve_size = 0;   
                    for(size_t j = 0; j < levels; j++)
                    {
                        std::cout << "Retrieve " << +level_num_bitplanes[j] << " (" << +(level_num_bitplanes[j] - prev_level_num_bitplanes[j]) << " more) bitplanes from level " << j << std::endl;
                        adios2::Engine one_level_values_reader_engine =
                            reader_io.Open(level_locations[j], adios2::Mode::Read);  
                        std::cout << "Retrieve from " << level_locations[j] << std::endl;
                        std::string varOneLevelValuesName = variableName+":Values:level-"+std::to_string(j);
                        auto varOneLevelValues = reader_io.InquireVariable<uint8_t>(varOneLevelValuesName);
                        varOneLevelValues.SetSelection({{offsets[j]}, {retrieve_sizes[j]}});
                        //std::vector<uint8_t> one_level_values(retrieve_sizes[j]);
                        uint8_t * buffer = (uint8_t *) malloc(retrieve_sizes[j]);
                        one_level_values_reader_engine.Get(varOneLevelValues, buffer, adios2::Mode::Sync);                    
                        //uint8_t * buffer = (uint8_t *) malloc(retrieve_sizes[j]);
                        //std::copy(one_level_values.begin()+offsets[j], one_level_values.begin()+offsets[j]+retrieve_sizes[j], buffer);
                        concated_level_components.push_back(buffer);
                        offsets[j] += retrieve_sizes[j];
                        total_retrieve_size += offsets[j];
                        one_level_values_reader_engine.Close();
                    }    
                    std::cout << "Total retrieve size = " << total_retrieve_size << std::endl;
                    for(size_t j = 0; j < level_num_bitplanes.size(); j++)
                    {
                        const uint8_t * pos = concated_level_components[j];
                        std::vector<const uint8_t*> interleaved_level;
                        for(size_t k = prev_level_num_bitplanes[j]; k < level_num_bitplanes[j]; k++){
                            interleaved_level.push_back(pos);
                            pos += level_sizes[j][k];
                        }
                        level_components.push_back(interleaved_level);
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
                    reconstructedData.clear();
                    reconstructedData = std::vector<T>(num_elements, 0);
                    auto level_elements = MDR::compute_level_elements(level_dims, target_level);
                    std::vector<uint32_t> dims_dummy(reconstruct_dimensions.size(), 0);
                    for(size_t j = 0; j <= target_level; j++)
                    {
                        compressor.decompress_level(level_components[j], level_sizes[j], prev_level_num_bitplanes[j], level_num_bitplanes[j] - prev_level_num_bitplanes[j]);
                        int level_exp = 0;
                        frexp(level_error_bounds[j], &level_exp);
                        auto level_decoded_data = encoder.progressive_decode(level_components[j], level_elements[j], level_exp, prev_level_num_bitplanes[j], level_num_bitplanes[j] - prev_level_num_bitplanes[j], j);
                        compressor.decompress_release();
                        const std::vector<uint32_t>& prev_dims = (j == 0) ? dims_dummy : level_dims[j-1];
                        interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[j], prev_dims, reconstructedData.data());
                        free(level_decoded_data);
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
                        exit(0);
                    }
                    for(size_t j = 0; j < concated_level_components.size(); j++)
                    {
                        free(concated_level_components[j]);
                    }
                    concated_level_components.clear();

                    err = clock_gettime(CLOCK_REALTIME, &end);
                    std::cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << std::endl;  
                    MGARD::print_statistics(rawVariableData.data(), reconstructedData.data(), rawVariableData.size());                                 
                }
            }          
        }
              
    }
    

    metadata_reader_engine.Close();
    
}