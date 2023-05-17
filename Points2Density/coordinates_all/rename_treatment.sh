#!/bin/bash                                                                     

d=$1

cp -- ~/Desktop/spatial_analysis/res_coord_series_2_scaled/$d/{A,B,C,D}*.txt $d/.

mv $d/A${d}-C.tif.txt $d/Pac0_5_R1-$d.txt
mv $d/A${d}-E.tif.txt $d/Pac0_5_R2-$d.txt
mv $d/A${d}-N.tif.txt $d/Pac0_5_R3-$d.txt
mv $d/A${d}-S.tif.txt $d/Pac0_5_R4-$d.txt
mv $d/A${d}-W.tif.txt $d/Pac0_5_R5-$d.txt

mv $d/B${d}-C.tif.txt $d/Pac0_05_R1-$d.txt
mv $d/B${d}-E.tif.txt $d/Pac0_05_R2-$d.txt
mv $d/B${d}-N.tif.txt $d/Pac0_05_R3-$d.txt
mv $d/B${d}-S.tif.txt $d/Pac0_05_R4-$d.txt
mv $d/B${d}-W.tif.txt $d/Pac0_05_R5-$d.txt

mv $d/C${d}-C.tif.txt $d/Pac0_005_R1-$d.txt
mv $d/C${d}-E.tif.txt $d/Pac0_005_R2-$d.txt
mv $d/C${d}-N.tif.txt $d/Pac0_005_R3-$d.txt
mv $d/C${d}-S.tif.txt $d/Pac0_005_R4-$d.txt
mv $d/C${d}-W.tif.txt $d/Pac0_005_R5-$d.txt

mv $d/D${d}-C.tif.txt $d/Pac0_0005_R1-$d.txt
mv $d/D${d}-E.tif.txt $d/Pac0_0005_R2-$d.txt
mv $d/D${d}-N.tif.txt $d/Pac0_0005_R3-$d.txt
mv $d/D${d}-S.tif.txt $d/Pac0_0005_R4-$d.txt
mv $d/D${d}-W.tif.txt $d/Pac0_0005_R5-$d.txt






                                                   
