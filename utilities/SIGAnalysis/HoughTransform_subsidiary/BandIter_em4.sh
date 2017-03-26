#!/bin/bash

for (( i=9; i<11; i++ ))
do
    root -b -l -q "HoughTransform_BandIter_em.C($i)"
done
