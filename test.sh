#!/bin/bash
make build
./build/debug/test/utest --gtest_filter="-test_large*"
RESULT=$?
[ $RESULT -ne 0 ] && exit 1
exit 0
