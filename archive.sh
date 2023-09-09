#!/bin/bash

set -ex 

make clean
make single-header gen-transpose-autogen format doc
touch T-SIMD_CODE.tar.gz # prevent tar from complaining about files changing while it is running
GZIP=-9 tar -vcz --exclude-ignore=.tarignore -f T-SIMD_CODE.tar.gz ../T-SIMD_CODE
