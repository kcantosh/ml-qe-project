#!/bin/bash

more pred0.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred0_clean.txt
more pred1.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred1_clean.txt
more pred2.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred2_clean.txt
more pred3.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred3_clean.txt
more pred4.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred4_clean.txt
more pred5.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred5_clean.txt
more pred6.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred6_clean.txt
more pred7.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred7_clean.txt
more pred8.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred8_clean.txt
more pred9.csv | awk 'BEGIN{FS=","}{if (NR>1){print $2,$3,$4,$5,$6} }' > pred9_clean.txt




