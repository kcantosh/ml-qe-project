more a_svm1000.txt | awk 'BEGIN {FS=":"} {for (i=1; i<130; i++) print $i}' | awk '{print $1}' | sed 's/\ //g' > a_data.txt
more c_svm1000.txt | awk 'BEGIN {FS=":"} {for (i=1; i<130; i++) print $i}' | awk '{print $1}' | sed 's/\ //g' > c_data.txt
more Ox_svm1000.txt | awk 'BEGIN {FS=":"} {for (i=1; i<130; i++) print $i}' | awk '{print $1}' | sed 's/\ //g' > Ox_data.txt
more Oy_svm1000.txt | awk 'BEGIN {FS=":"} {for (i=1; i<130; i++) print $i}' | awk '{print $1}' | sed 's/\ //g' > Oy_data.txt
more Oz_svm1000.txt | awk 'BEGIN {FS=":"} {for (i=1; i<130; i++) print $i}' | awk '{print $1}' | sed 's/\ //g' > Oz_data.txt


