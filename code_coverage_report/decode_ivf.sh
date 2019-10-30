#!/bin/bash

# find all ivf files to decoder 
for file_list in $(ls .)
do
    if [ ${file_list##*.} = "ivf" ]; then
        # setup for output YUV file
        output_yuv=${file_list%.*}.yuv
		./SvtAv1DecApp -i ${file_list} -o ${output_yuv}
        # delete the output YUV file to save disk space
        rm -f ${output_yuv}
    fi
done
