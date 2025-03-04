#!/bin/bash

./model/bin/w3_clean -m
./model/bin/w3_setup model -c intel -s $1
./model/bin/w3_automake

