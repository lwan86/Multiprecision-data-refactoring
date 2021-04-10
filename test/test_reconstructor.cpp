#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <iomanip>
#include <cmath>
#include <bitset>
#include "utils.hpp"
#include "Reconstructor/Reconstructor.hpp"
#include <mpi.h>
// #include "evaluate.hpp"

using namespace std;

template <class T, class Reconstructor>
void evaluate(const vector<T>& data, const vector<double>& tolerance, Reconstructor reconstructor){
    struct timespec start, end;
    int err = 0;
    auto dims = reconstructor.get_dimensions();
    double time = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for(int i=0; i<tolerance.size(); i++){
        // cout << "Start reconstruction" << endl;
        // err = clock_gettime(CLOCK_REALTIME, &start);
        auto reconstructed_data = reconstructor.progressive_reconstruct(tolerance[i]);
        // err = clock_gettime(CLOCK_REALTIME, &end);
        // cout << "Reconstruct time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_nsec - start.tv_nsec)/(double)1000000000 << "s" << endl;
        // TODO: add full resolution check
        if(rank == 0) MGARD::print_statistics(data.data(), reconstructed_data, data.size());
        // COMP_UTILS::evaluate_gradients(data.data(), reconstructed_data, dims[0], dims[1], dims[2]);
        // COMP_UTILS::evaluate_average(data.data(), reconstructed_data, dims[0], dims[1], dims[2], 5);
    }
}

template <class T, class Decomposer, class Interleaver, class Encoder, class Compressor, class ErrorEstimator, class SizeInterpreter, class Retriever>
void test(string filename, vector<double>& tolerance, Decomposer decomposer, Interleaver interleaver, Encoder encoder, Compressor compressor, ErrorEstimator estimator, SizeInterpreter interpreter, Retriever retriever){
    auto reconstructor = MDR::ComposedReconstructor<T, Decomposer, Interleaver, Encoder, Compressor, SizeInterpreter, ErrorEstimator, Retriever>(decomposer, interleaver, encoder, compressor, interpreter, retriever);
    // cout << "loading metadata" << endl;
    double time = 0;
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        time = - MPI_Wtime();
    }
    reconstructor.load_metadata();
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        time = + MPI_Wtime();
        cout << "metadata reading time = " << time << endl;
    }    
    size_t num_elements = 0;
    auto data = MGARD::readfile<T>(filename.c_str(), num_elements);
    // change tolerance to PSNR
    {
        T max_val = data[0];
        T min_val = data[0];
        for(int i=0; i<data.size(); i++){
            if(max_val < data[i]) max_val = data[i];
            if(min_val > data[i]) min_val = data[i];
        }
        T value_range = max_val - min_val;
        for(int i=0; i<tolerance.size(); i++){
            tolerance[i] = (value_range / pow(10, tolerance[i]/20))*(value_range / pow(10, tolerance[i]/20)) * data.size() * size;        
        }
    }
    evaluate(data, tolerance, reconstructor);
}

int main(int argc, char ** argv){
    MPI_Init(&argc, &argv);
    int argv_id = 1;
    string filename = string(argv[argv_id ++]);
    int error_mode = atoi(argv[argv_id++]);
    int num_tolerance = atoi(argv[argv_id ++]);
    vector<double> tolerance(num_tolerance, 0);
    for(int i=0; i<num_tolerance; i++){
        tolerance[i] = atof(argv[argv_id ++]);    
    }
    double s = atof(argv[argv_id ++]);

    string metadata_file = "refactored_data/metadata.bin";
    int num_levels = 0;
    int num_dims = 0;
    {
        // metadata interpreter, otherwise information needs to be provided
        // size_t num_bytes = 0;
        // auto metadata = MGARD::readfile<uint8_t>(metadata_file.c_str(), num_bytes);
        // assert(num_bytes > num_dims * sizeof(uint32_t) + 2);
        // num_dims = metadata[0];
        // num_levels = metadata[num_dims * sizeof(uint32_t) + 1];
        // cout << "number of dimension = " << num_dims << ", number of levels = " << num_levels << endl;

        // fix num_dims and num_levels for this example
        num_dims = 3;
        num_levels = 4;
    }
    vector<string> files;
    for(int i=0; i<num_levels; i++){
        string filename = "refactored_data/level_" + to_string(i) + ".bin";
        files.push_back(filename);
    }

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
    // auto retriever = MDR::ConcatLevelFileRetriever(metadata_file, files);
    auto retriever = MDR::HPSSFileRetriever(metadata_file, files);
    switch(error_mode){
        case 1:{
            auto estimator = MDR::SNormErrorEstimator<T>(num_dims, num_levels - 1, s);
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            // auto interpreter = MDR::InorderSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            // auto interpreter = MDR::RoundRobinSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            auto interpreter = MDR::NegaBinaryGreedyBasedSizeInterpreter<MDR::SNormErrorEstimator<T>>(estimator);
            // auto estimator = MDR::L2ErrorEstimator_HB<T>(num_dims, num_levels - 1);
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::L2ErrorEstimator_HB<T>>(estimator);
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);            
            break;
        }
        default:{
            auto estimator = MDR::MaxErrorEstimatorOB<T>(num_dims);
            auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorOB<T>>(estimator);
            // auto estimator = MDR::MaxErrorEstimatorHB<T>();
            // auto interpreter = MDR::SignExcludeGreedyBasedSizeInterpreter<MDR::MaxErrorEstimatorHB<T>>(estimator);
            test<T>(filename, tolerance, decomposer, interleaver, encoder, compressor, estimator, interpreter, retriever);
        }
    }    
    MPI_Finalize();
    return 0;
}
