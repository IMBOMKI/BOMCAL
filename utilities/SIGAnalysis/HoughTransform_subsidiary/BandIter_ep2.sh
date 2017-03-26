#!/bin/bash

for (( i=7; i<11; i++ ))
do
    root -b -l -q "HoughTransform_BandIter_ep.C($i)"
done
