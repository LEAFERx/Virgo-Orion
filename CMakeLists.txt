cmake_minimum_required(VERSION 3.16)

project(superorion)

set(CMAKE_C_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)

include_directories(lib)
include_directories(include)

add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -O3")
#we enable optimization by default, you can remove -O3 for debug purpose
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -lstdc++ -Wall -march=native -faligned-new -O3")

link_directories(lib)

aux_source_directory(src/common COMMON_SRC)

# Orion

add_executable(linear_code
    ${COMMON_SRC}
    src/orion/test/encode_test.cpp
    src/orion/linear_code/linear_code_encode.cpp
    src/orion/linear_code/expander.cpp
    src/orion/infrastructure/RS_polynomial.cpp
    src/orion/infrastructure/merkle_tree.cpp)

add_executable(linearPC
    ${COMMON_SRC}
    src/orion/test/univariate_PC_test.cpp
    src/orion/linear_code/linear_code_encode.cpp
    src/orion/linear_code/expander.cpp
    src/orion/VPD/linearPC.cpp
    src/orion/infrastructure/RS_polynomial.cpp
    src/orion/infrastructure/merkle_tree.cpp
    src/orion/VPD/fri.cpp
    src/orion/VPD/vpd_prover.cpp
    src/orion/VPD/vpd_verifier.cpp
    src/orion/linear_gkr/verifier.cpp
    src/orion/linear_gkr/prover.cpp
    src/orion/poly_commitment/poly_commit.cpp)

add_executable(linearPC_multi
    ${COMMON_SRC}
    src/orion/test/multivariate_PC_test.cpp
    src/orion/linear_code/linear_code_encode.cpp
    src/orion/linear_code/expander.cpp
    src/orion/VPD/linearPC.cpp
    src/orion/infrastructure/RS_polynomial.cpp
    src/orion/infrastructure/merkle_tree.cpp
    src/orion/VPD/fri.cpp
    src/orion/VPD/vpd_prover.cpp
    src/orion/VPD/vpd_verifier.cpp
    src/orion/linear_gkr/verifier.cpp
    src/orion/linear_gkr/prover.cpp
    src/orion/poly_commitment/poly_commit.cpp)

target_link_libraries(linear_code -Llib -lflo-shani -lXKCP -lgmp)
target_link_libraries(linearPC -Llib -lflo-shani -lXKCP -lgmp)
target_link_libraries(linearPC_multi -Llib -lflo-shani -lXKCP -lgmp)

aux_source_directory(src/orion/infrastructure ORION_INFRA_SRC)
aux_source_directory(src/orion/linear_code ORION_LCODE_SRC)
aux_source_directory(src/orion/linear_gkr ORION_LGKR_SRC)
aux_source_directory(src/orion/poly_commitment ORION_PC_SRC)
aux_source_directory(src/orion/VPD ORION_VPD_SRC)
add_library(orion ${COMMON_SRC} ${ORION_INFRA_SRC} ${ORION_LCODE_SRC} ${ORION_LGKR_SRC} ${ORION_PC_SRC} ${ORION_VPD_SRC})
target_link_libraries(orion -Llib -lflo-shani -lXKCP -lgmp)

# Frontend

add_executable(fft_gkr
    ${COMMON_SRC}
    src/frontend/fft_circuit_GKR.cpp)

add_executable(zk_proof
    ${COMMON_SRC}
    src/frontend/zk_verifier.cpp
    src/frontend/zk_prover.cpp
    src/frontend/main_zk.cpp)

target_link_libraries(fft_gkr -lgmp)

target_link_libraries(zk_proof -Llib -lflo-shani -lXKCP -lgmp orion)