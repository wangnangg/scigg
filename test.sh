#!/bin/bash
make build
./build/debug/test/utest
RESULT=$?
[ $RESULT -ne 0 ] && exit 1
exit 1
