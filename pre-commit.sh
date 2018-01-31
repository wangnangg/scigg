#!/bin/bash
./format.sh
git add .
./test.sh debug "-test_large*"
