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

int main(int argc, char ** argv)
{
    int argv_id = 1;
    std::string refactoredFileName = std::string(argv[argv_id ++]);
    int error_mode = atoi(argv[argv_id++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    std::vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);    
    }
    double s = atof(argv[argv_id ++]);    
    std::string originalFileName = std::string(argv[argv_id ++]);

    adios2::ADIOS adios;
    adios2::IO reader_io = adios.DeclareIO("ReaderIO");
    adios2::Engine reader_engine =
           reader_io.Open(refactoredFileName, adios2::Mode::Read);  
    auto varDimensions = reader_io.InquireVariable<uint32_t>("Dimensions");
    std::vector<uint32_t> dimensions(varDimensions.Shape()[0]);
    reader_engine.Get(varDimensions, dimensions.data(), adios2::Mode::Sync);
    std::cout << "dimensions: ";
    for (size_t i = 0; i < dimensions.size(); i++)
    {
        std::cout << dimensions[i] << " ";
    }
    std::cout << std::endl;
    uint32_t levels;
    auto varNumLevels = reader_io.InquireVariable<uint32_t>("Levels");
    reader_engine.Get(varNumLevels, levels, adios2::Mode::Sync);
    std::cout << "# of levels: " << levels << std::endl;

    using T = float;
    using T_stream = uint32_t;
    std::vector<T> level_error_bounds(levels);
    std::vector<std::vector<double>> level_squared_errors(levels);
    std::vector<std::vector<uint32_t>> level_sizes(levels);
    std::vector<uint8_t> level_num_bitplanes;
    std::vector<uint8_t> stopping_indices(levels);
    level_num_bitplanes = std::vector<uint8_t>(levels, 0);
    for (size_t i = 0; i < levels; i++)
    {
        std::string varOneLevelSizesName = "/Sizes/level-"+std::to_string(i);
        auto varOneLevelSizes = reader_io.InquireVariable<uint32_t>(varOneLevelSizesName);
        std::vector<uint32_t> one_level_sizes(varOneLevelSizes.Shape()[0]);
        reader_engine.Get(varOneLevelSizes, one_level_sizes.data(), adios2::Mode::Sync);
        level_sizes[i] = one_level_sizes;

        std::string varOneLevelSquaredErrorsName = "/SquaredErrors/level-"+std::to_string(i);
        auto varOneLevelSquaredErrors = reader_io.InquireVariable<double>(varOneLevelSquaredErrorsName);
        std::vector<double> one_level_squared_errors(varOneLevelSquaredErrors.Shape()[0]);
        reader_engine.Get(varOneLevelSquaredErrors, one_level_squared_errors.data(), adios2::Mode::Sync);
        level_squared_errors[i] = one_level_squared_errors;

        std::string varOneLevelErrorBoundsName = "/ErrorBounds/level-"+std::to_string(i);
        auto varOneLevelErrorBounds = reader_io.InquireVariable<T>(varOneLevelErrorBoundsName);
        T one_level_error_bounds;
        reader_engine.Get(varOneLevelErrorBounds, one_level_error_bounds, adios2::Mode::Sync);
        level_error_bounds[i] = one_level_error_bounds;

        std::string varOneLevelStopIndicesName = "/StopIndices/level-"+std::to_string(i);
        auto varOneLevelStopIndices = reader_io.InquireVariable<uint8_t>(varOneLevelStopIndicesName);
        uint8_t one_level_stop_index;
        reader_engine.Get(varOneLevelStopIndices, one_level_stop_index, adios2::Mode::Sync);
        stopping_indices[i] = one_level_stop_index;
    }
    // for (size_t i = 0; i < levels; i++)
    // {
    //     std::cout << "/Sizes/level-" << i << " has "<< level_sizes[i].size() << " elements: ";
    //     for (size_t j = 0; j < level_sizes[i].size(); j++)
    //     {
    //         std::cout << level_sizes[i][j] << " ";
    //     }
    //     std::cout << std::endl;  
    //     std::cout << "/SquaredErrors/level-" << i << " has "<< level_squared_errors[i].size() << " elements: ";
    //     for (size_t j = 0; j < level_squared_errors[i].size(); j++)
    //     {
    //         std::cout << level_squared_errors[i][j] << " ";
    //     }
    //     std::cout << std::endl;  
    //     std::cout << "/ErrorBounds/level-" << i << ": ";
    //     std::cout << level_error_bounds[i] << std::endl;   
    // }    

    size_t original_num_elements = 0;
    auto originalData = MGARD::readfile<T>(originalFileName.c_str(), original_num_elements);    

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
    struct timespec start, end;
    int err = 0;
    switch(error_mode)
    {
        case 1:
        {
            auto estimator = MDR::SNormErrorEstimator<T>(dimensions.size(), levels - 1, s);
            auto interpreter = MDR::NegaBinaryGreedyBasedSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);

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
                for (size_t j = 0; j < level_sizes.size(); j++)
                {
                    std::cout << j << ": ";
                    for (size_t k = 0; k < level_sizes[j].size(); k++)
                    {
                        std::cout << level_sizes[j][k] << " ";
                    }
                    std::cout << std::endl;
                }
                for (size_t j = 0; j < level_num_bitplanes.size(); j++)
                {
                    std::cout << +level_num_bitplanes[j] << " ";
                }
                std::cout << std::endl;

                auto retrieve_sizes = interpreter.interpret_retrieve_size(level_sizes, level_errors, tolerance[i], level_num_bitplanes);
                std::vector<std::vector<const uint8_t*>> level_components;
                std::vector<uint8_t*> concated_level_components;
                uint32_t total_retrieve_size = 0; 
                for(size_t j = 0; j < levels; j++)
                {
                    std::cout << "Retrieve " << +level_num_bitplanes[j] << " (" << +(level_num_bitplanes[j] - prev_level_num_bitplanes[j]) << " more) bitplanes from level " << j << std::endl;
                    std::string varOneLevelValueName = "/Value/level-"+std::to_string(j);
                    auto varOneLevelValue = reader_io.InquireVariable<uint8_t>(varOneLevelValueName);
                    std::vector<uint8_t> one_level_value(varOneLevelValue.Shape()[0]);
                    reader_engine.Get(varOneLevelValue, one_level_value.data(), adios2::Mode::Sync);                    
                    uint8_t * buffer = (uint8_t *) malloc(retrieve_sizes[j]);
                    std::copy(one_level_value.begin()+offsets[j], one_level_value.begin()+offsets[j]+retrieve_sizes[j], buffer);
                    concated_level_components.push_back(buffer);
                    offsets[j] += retrieve_sizes[j];
                    total_retrieve_size += offsets[j];
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
                    //std::cout << "decompress_level" << std::endl;
                    int level_exp = 0;
                    frexp(level_error_bounds[j], &level_exp);
                    //std::cout << "frexp" << std::endl;
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
                    //std::cout << "progressive_decode" << std::endl;
                    compressor.decompress_release();
                    //std::cout << "decompress_release" << std::endl;
                    const std::vector<uint32_t>& prev_dims = (j == 0) ? dims_dummy : level_dims[j-1];
                    interleaver.reposition(level_decoded_data, reconstruct_dimensions, level_dims[j], prev_dims, reconstructedData.data());
                    //std::cout << "reposition" << std::endl;
                    free(level_decoded_data);
                    //std::cout << "free" << std::endl;
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
                MGARD::print_statistics(originalData.data(), reconstructedData.data(), originalData.size());               
            }
                     
        }
    }
    reader_engine.Close();
    return 0;
}