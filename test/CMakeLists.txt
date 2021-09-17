set(EVA_INCLUDES "../external/compression_utils/include")

add_executable (test_encoder test_encoder.cpp)
target_include_directories(test_encoder PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_encoder ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_decomposer test_decomposer.cpp)
target_include_directories(test_decomposer PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_decomposer ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_interleaver test_interleaver.cpp)
target_include_directories(test_interleaver PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_interleaver ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_error_collector test_error_collector.cpp)
target_include_directories(test_error_collector PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_error_collector ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_refactor test_refactor.cpp)
target_include_directories(test_refactor PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_refactor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})

add_executable (test_reconstructor test_reconstructor.cpp)
target_include_directories(test_reconstructor PRIVATE ${EVA_INCLUDES} ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES})
target_link_libraries(test_reconstructor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB})


#add_executable (test_refactor_adios2 test_refactor_adios2.cpp)
#target_include_directories(test_refactor_adios2 PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${ADIOS2_INCLUDES})
#target_link_libraries(test_refactor_adios2 ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${ADIOS2_LIB})

#add_executable (test_reconstructor_adios2 test_reconstructor_adios2.cpp)
#target_include_directories(test_reconstructor_adios2 PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${ADIOS2_INCLUDES})
#target_link_libraries(test_reconstructor_adios2 ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${ADIOS2_LIB})

add_executable (sirius2_data_refactor sirius2_data_refactor.cpp)
target_include_directories(sirius2_data_refactor PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${ADIOS2_INCLUDES})
target_link_libraries(sirius2_data_refactor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${ADIOS2_LIB})

add_executable (sirius2_data_reconstructor sirius2_data_reconstructor.cpp)
target_include_directories(sirius2_data_reconstructor PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${ADIOS2_INCLUDES})
target_link_libraries(sirius2_data_reconstructor ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${ADIOS2_LIB})

add_executable (sirius2_data_refactor_ec sirius2_data_refactor_ec.cpp)
target_include_directories(sirius2_data_refactor_ec PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${ADIOS2_INCLUDES} ${EC_INCLUDES})
target_link_libraries(sirius2_data_refactor_ec ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${ADIOS2_LIB} ${EC_LIB})

add_executable (sirius2_data_reconstructor_ec sirius2_data_reconstructor_ec.cpp)
target_include_directories(sirius2_data_reconstructor_ec PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${ADIOS2_INCLUDES} ${EC_INCLUDES})
target_link_libraries(sirius2_data_reconstructor_ec ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${ADIOS2_LIB} ${EC_LIB})


add_executable (test_mgard_ec_mpi test_mgard_ec_mpi.cpp)
message(${ADIOS2MPI_LIB})
message(${ADIOS2MPI_INCLUDES})
target_include_directories(test_mgard_ec_mpi PRIVATE ${MGARDx_INCLUDES} ${SZ3_INCLUDES} ${ZSTD_INCLUDES} ${ADIOS2_INCLUDES} ${ADIOS2MPI_INCLUDES} ${EC_INCLUDES} ${MPI_INCLUDES} ${ROCKSDB_INCLUDES})
target_link_libraries(test_mgard_ec_mpi ${PROJECT_NAME} ${SZ3_LIB} ${ZSTD_LIB} ${ADIOS2_LIB} ${ADIOS2MPI_LIB} ${EC_LIB} ${MPI_LIB} ${ROCKSDB_LIB})
target_compile_definitions(test_mgard_ec_mpi PRIVATE ADIOS2_USE_MPI)