#!/bin/bash

for f in *.JOB
do
    java -jar ../dist/SRA_DNLD_MGR.jar "$f"
    rc=$?; if [[ $rc -ne 0 ]]; then echo "error job $f"; exit $rc; fi
done
