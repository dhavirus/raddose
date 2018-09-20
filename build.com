#!/bin/bash

g77 -O -w -fno-second-underscore -fno-globals -fno-automatic -o raddose name2z.f z2name.f mucal.f upcase.f raddose.f -L/programs/CCP4/ccp4_6.1.13/ccp4-6.1.13/lib/libccp4f.a
