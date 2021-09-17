#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>

#include "../include/Decomposer/Decomposer.hpp"
#include "../include/Interleaver/Interleaver.hpp"
#include "../include/BitplaneEncoder/BitplaneEncoder.hpp"
#include "../include/ErrorCollector/ErrorCollector.hpp"
#include "../include/LosslessCompressor/LevelCompressor.hpp"
#include "../include/RefactorUtils.hpp"

#include <adios2.h>

int main(int argc, char *argv[])
{
    std::string inputFileName;
    size_t storageTiers = 0;
    std::vector<std::string> storageTiersPaths;
    size_t target_level = 0;
    size_t num_bitplanes = 0;

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
                std::stringstream ssStorageTiers(argv[i+1]);
                ssStorageTiers >> storageTiers;
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
                        storageTiersPaths.push_back(argv[j]);
                    }  
                }
                else
                {
                    std::cerr << "--tiers option requires [# of tiers] arguments." << std::endl;
                    return 1;
                } 
            }
            else
            {
                std::cerr << "--[# of tiers] needs to be greater than zero." << std::endl;
                return 1;
            }
        }
        else if (arg == "-l" || arg == "--levels")
        {
            if (i+1 < argc)
            {
                std::stringstream ssLevels(argv[i+1]);
                ssLevels >> target_level;
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
                std::stringstream ssPlanes(argv[i+1]);
                ssPlanes >> num_bitplanes;
            }
            else
            {
                std::cerr << "--planes option requires one argument." << std::endl;
                return 1;
            }            
        }                
    } 

    adios2::ADIOS adios;
    adios2::IO writer_io = adios.DeclareIO("WriterIO");
    std::vector<adios2::Engine> data_writer_engines(storageTiersPaths.size());    
    for (size_t i = 0; i < storageTiersPaths.size(); i++)
    {
        std::cout << "tier " << i << ": " << storageTiersPaths[i] << std::endl;
        std::string fullDataPath = storageTiersPaths[i];
        std::string refactoredDataFileName = "refactored.data.bp";
        if (!fullDataPath.empty() && fullDataPath.back() != '/')
        {
            fullDataPath += '/';
        }
        fullDataPath += refactoredDataFileName;
        //std::cout << fullDataPath << std::endl;
        adios2::Engine data_writer_engine =
           writer_io.Open(fullDataPath, adios2::Mode::Write); 
        data_writer_engines[i] = data_writer_engine;
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
        if (variableType == "float")
        {
            auto variable = reader_io.InquireVariable<float>(variableName);
            size_t spaceDimensions = variable.Shape().size();
            if (spaceDimensions == 1)
            {
                continue;
            }
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
            // auto encoder = MDR::GroupedBPEncoder<T, T_stream>();
            auto encoder = MDR::PerBitBPEncoder<T, T_stream>();
            auto compressor = MDR::DefaultLevelCompressor();
            // auto compressor = MDR::NullLevelCompressor();
            auto collector = MDR::SquaredErrorCollector<T>();
            //std::cout << num_bitplanes << std::endl;

            std::vector<uint32_t> dimensions(variable.Shape().size());
            std::vector<T> level_error_bounds;
            std::vector<std::vector<uint8_t*>> level_components;
            std::vector<std::vector<uint32_t>> level_sizes;
            std::vector<std::vector<double>> level_squared_errors;
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
                auto collected_error = s_collector.collect_level_error(buffer, level_elements[i], num_bitplanes, level_max_error);
                // std::cout << collected_error.size() << std::endl;
                level_squared_errors.push_back(collected_error);
                // encode level data
                int level_exp = 0;
                frexp(level_max_error, &level_exp);
                std::vector<uint32_t> stream_sizes;
                auto streams = encoder.encode(buffer, level_elements[i], level_exp, num_bitplanes, stream_sizes);
                free(buffer);
                // lossless compression
                compressor.compress_level(streams, stream_sizes);
                // record encoded level data and size
                level_components.push_back(streams);
                level_sizes.push_back(stream_sizes);
            }

            for (size_t i = 0; i < level_components.size(); i++)
            {
                uint32_t concated_level_size = 0;
                for(size_t j = 0; j < level_components[i].size(); j++)
                {
                    concated_level_size += level_sizes[i][j];
                }
                uint8_t * concated_level_data = (uint8_t *) malloc(concated_level_size);
                uint8_t * concated_level_data_pos = concated_level_data;
                for(size_t j = 0; j < level_components[i].size(); j++)
                {
                    memcpy(concated_level_data_pos, level_components[i][j], level_sizes[i][j]);
                    concated_level_data_pos += level_sizes[i][j];
                } 
                std::string varLevelValuesName = variableName+":Values:level-"+std::to_string(i); 
                int assignedTier = static_cast<int>(static_cast<double>(i)*(static_cast<double>(static_cast<double>(storageTiers)/level_components.size())));
                //std::cout << "level " << i << ", " << "tier " << assignedTier << std::endl;
                adios2::Variable<uint8_t> varLevelValues = writer_io.DefineVariable<uint8_t>(varLevelValuesName, {concated_level_size}, {0}, {concated_level_size});
                data_writer_engines[assignedTier].Put(varLevelValues, concated_level_data, adios2::Mode::Sync);
                std::string varLevelSizesName = variableName+":Sizes:level-"+std::to_string(i);
                adios2::Variable<uint32_t> varLevelSizes = writer_io.DefineVariable<uint32_t>(varLevelSizesName, {level_sizes[i].size()}, {0}, {level_sizes[i].size()});
                metadata_writer_engine.Put(varLevelSizes, level_sizes[i].data(), adios2::Mode::Sync);   
                std::string varLevelLocationsName = variableName+":Locations:level-"+std::to_string(i);
                adios2::Variable<std::string> varLevelLocations = writer_io.DefineVariable<std::string>(varLevelLocationsName);
                std::string fullDataPath = storageTiersPaths[assignedTier];
                std::string refactoredDataFileName = "refactored.data.bp";
                if (!fullDataPath.empty() && fullDataPath.back() != '/')
                {
                    fullDataPath += '/';
                }
                fullDataPath += refactoredDataFileName;
                metadata_writer_engine.Put(varLevelLocations, fullDataPath, adios2::Mode::Sync);         
            }
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

            for (size_t i = 0; i < level_error_bounds.size(); i++)
            {
                std::string varLevelErrorBoundsName = variableName+":ErrorBounds:level-"+std::to_string(i);
                adios2::Variable<T> varLevelErrorBounds = writer_io.DefineVariable<T>(varLevelErrorBoundsName);
                metadata_writer_engine.Put(varLevelErrorBounds, level_error_bounds[i], adios2::Mode::Sync);
                std::string varLevelSquaredErrorsName = variableName+":SquaredErrors:level-"+std::to_string(i);
                adios2::Variable<double> varLevelSquaredErrors = writer_io.DefineVariable<double>(varLevelSquaredErrorsName, {level_squared_errors[i].size()}, {0}, {level_squared_errors[i].size()});
                metadata_writer_engine.Put(varLevelSquaredErrors, level_squared_errors[i].data(), adios2::Mode::Sync);
            }
        } 
    }    
    for (size_t i = 0; i < data_writer_engines.size(); i++)
    {
        data_writer_engines[i].Close();
    }
    metadata_writer_engine.Close();
    
    reader_engine.Close();
}