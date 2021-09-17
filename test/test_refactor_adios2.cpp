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

using namespace std;


int main(int argc, char ** argv){

    int argv_id = 1;
    std::string inputFileName = std::string(argv[argv_id ++]);
    int target_level = atoi(argv[argv_id ++]);
    int num_bitplanes = atoi(argv[argv_id ++]);
    int num_dims = atoi(argv[argv_id ++]);
    std::vector<uint32_t> dims(num_dims, 0);
    for(size_t i = 0; i < num_dims; i++)
    {
        dims[i] = atoi(argv[argv_id ++]);
    }
    std::string outputFileName = std::string(argv[argv_id ++]);

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
    auto collector = MDR::SquaredErrorCollector<T>();

    std::vector<uint32_t> dimensions;
    std::vector<T> level_error_bounds;
    std::vector<std::vector<uint8_t*>> level_components;
    std::vector<std::vector<uint32_t>> level_sizes;
    std::vector<std::vector<double>> level_squared_errors;
    std::vector<uint8_t> stopping_indices;

    adios2::ADIOS adios;
    adios2::IO writer_io = adios.DeclareIO("WriterIO");
    adios2::Engine writer_engine =
           writer_io.Open(outputFileName, adios2::Mode::Write);   
    
    size_t n_elements = 0;
    auto inputData = MGARD::readfile<T>(inputFileName.c_str(), n_elements);
    T const * data_ = inputData.data();
    struct timespec start, end;
    int err = 0;
    cout << "Start refactoring" << endl;
    err = clock_gettime(CLOCK_REALTIME, &start);

    dimensions = dims;
    uint32_t num_elements = 1;
    for(const auto& dim:dimensions)
    {
        num_elements *= dim;
    }
    std::vector<T> data;
    data = std::vector<T>(data_, data_ + num_elements);

    uint8_t max_level = log2(*min_element(dimensions.begin(), dimensions.end())) - 1;
    if(target_level > max_level)
    {
        std::cerr << "Target level is higher than " << max_level << std::endl;
        return 1;
    }
    // decompose data hierarchically
    decomposer.decompose(data.data(), dimensions, target_level);

    auto level_dims = MDR::compute_level_dims(dimensions, target_level);
    auto level_elements = MDR::compute_level_elements(level_dims, target_level);
    std::vector<uint32_t> dims_dummy(dimensions.size(), 0);
    //MDR::SquaredErrorCollector<T> s_collector = MDR::SquaredErrorCollector<T>();
    for(size_t i = 0; i <= target_level; i++)
    {
        const std::vector<uint32_t>& prev_dims = (i == 0) ? dims_dummy : level_dims[i - 1];
        T * buffer = (T *) malloc(level_elements[i] * sizeof(T));
        // extract level i component
        interleaver.interleave(data.data(), dimensions, level_dims[i], prev_dims, reinterpret_cast<T*>(buffer));
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

    std::string varNumLevelsName = "Levels";
    adios2::Variable<uint32_t> varNumLevels = writer_io.DefineVariable<uint32_t>(varNumLevelsName);
    uint32_t numLevels = level_components.size();
    writer_engine.Put(varNumLevels, numLevels, adios2::Mode::Sync);
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
        std::string varLevelName = "/Value/level-"+std::to_string(i);
        adios2::Variable<uint8_t> varLevel = writer_io.DefineVariable<uint8_t>(varLevelName, {concated_level_size}, {0}, {concated_level_size});
        writer_engine.Put(varLevel, concated_level_data, adios2::Mode::Sync);
        std::string varLevelSizesName = "/Sizes/level-"+std::to_string(i);
        adios2::Variable<uint32_t> varLevelSizes = writer_io.DefineVariable<uint32_t>(varLevelSizesName, {level_sizes[i].size()}, {0}, {level_sizes[i].size()});
        writer_engine.Put(varLevelSizes, level_sizes[i].data(), adios2::Mode::Sync);
    }

    // uint32_t metadata_size = sizeof(uint8_t) + get_size(dimensions) // dimensions
    //                         + sizeof(uint8_t) + get_size(level_error_bounds) + get_size(level_squared_errors) + get_size(level_sizes); // level information
    // uint8_t * metadata = (uint8_t *) malloc(metadata_size);
    // uint8_t * metadata_pos = metadata;
    // *(metadata_pos ++) = (uint8_t) dimensions.size();
    // MDR::serialize(dimensions, metadata_pos);
    // *(metadata_pos ++) = (uint8_t) level_error_bounds.size();
    // MDR::serialize(level_error_bounds, metadata_pos);
    // MDR::serialize(level_squared_errors, metadata_pos);
    // MDR::serialize(level_sizes, metadata_pos);
    std::string varDimensionName = "Dimensions";
    adios2::Variable<uint32_t> varDimension = writer_io.DefineVariable<uint32_t>(varDimensionName, {dimensions.size()}, {0}, {dimensions.size()});
    writer_engine.Put(varDimension, dimensions.data(), adios2::Mode::Sync);    

    for (size_t i = 0; i < level_error_bounds.size(); i++)
    {
        std::string varLevelErrorBoundName = "/ErrorBounds/level-"+std::to_string(i);
        adios2::Variable<T> varLevelErrorBound = writer_io.DefineVariable<T>(varLevelErrorBoundName);
        writer_engine.Put(varLevelErrorBound, level_error_bounds[i], adios2::Mode::Sync);
        std::string varLevelSquaredErrorName = "/SquaredErrors/level-"+std::to_string(i);
        adios2::Variable<double> varLevelSquaredError = writer_io.DefineVariable<double>(varLevelSquaredErrorName, {level_squared_errors[i].size()}, {0}, {level_squared_errors[i].size()});
        writer_engine.Put(varLevelSquaredError, level_squared_errors[i].data(), adios2::Mode::Sync);
        std::string varLevelStopIndicesName = "/StopIndices/level-"+std::to_string(i);
        adios2::Variable<uint8_t> varLevelStopIndices = writer_io.DefineVariable<uint8_t>(varLevelStopIndicesName);
        writer_engine.Put(varLevelStopIndices, stopping_indices[i], adios2::Mode::Sync);
    }

    writer_engine.Close();

    err = clock_gettime(CLOCK_REALTIME, &end);
    cout << "Refactor time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;

    return 0;
}