#!/bin/bash

for (( i=7; i<9; i++ ))
do
    root -b -l -q "HoughTransform_BandIter_em.C($i)"
done
