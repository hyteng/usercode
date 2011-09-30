#!/bin/bash
export FileList=`cat $1`
touch BatchValidation.bat
for file in ${FileList}
do
    sed -e "s|-file-|${file}|g" Validation.lsf > Validation_${file}.lsf
    chmod u+x Validation_${file}.lsf
    echo "bsub -q 1nh Validation_${file}.lsf" >> BatchValidation.bat
done
chmod u+x BatchValidation.bat
