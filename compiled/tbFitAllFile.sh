#!/bin/bash
# Batch-run tbFit.cc over all 11 PPM points (postMacroAllFile0..10.root).
# tbFit.cc always writes tbFitPPM10.00.root and canFitLevel*.pdf regardless of
# theFileNumber (dopant is hardcoded), so each run's output is renamed with
# the file index right after it completes to avoid overwriting the next run.

theFitChannel=-2

for i in 0 1 2 3 4 5 6 7 8 9 10
do
    root -b -q "tbFit.cc(${theFitChannel}, ${i})"

    for f in tbFitPPM10.00.root tbFitSimPPM10.00.root canFitLevel0.pdf canFitLevel1.pdf canFitLevel2.pdf canFitTrig.pdf
    do
        [ -f "$f" ] && mv "$f" "${f%.*}_file${i}.${f##*.}"
    done
done
