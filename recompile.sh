#!/bin/bash

autoreconf && ./configure && make clean && make -j

#cd src
#make check -j 
