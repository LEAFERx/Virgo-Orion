#!/bin/bash
cmake ../..
make -C ../.. zk_proof -j4
make -C ../.. fft_gkr -j4
mv ../../zk_proof .
mv ../../fft_gkr .