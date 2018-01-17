#!/bin/bash
config=$1
make utest config=${config}
./build/${config}/utest --gtest_filter="-test_large*"
RESULT=$?
[ $RESULT -ne 0 ] && exit 1
exit 0
