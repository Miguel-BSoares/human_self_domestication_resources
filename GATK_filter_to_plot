##script to extract information from different vcf files

#!/bin/bash

list="population_snps"

for file in $list
do
    
    zegrep -v "^#" ${file}.vcf.gz | \
    cut -f 8 | \
    sed 's/^.*;DP=\([0-9]*\);.*$/\1/' > ${file}_DP.txt
    
    
    zegrep -v "^#" ${file}.vcf.gz | \
    cut -f 8 | \
    sed 's/^.*;QD=\([0-9]*\.[0-9]*\);.*$/\1/' > ${file}_QD.txt
    
    
    zegrep -v "^#" ${file}.vcf.gz | \
    cut -f 8 | \
    sed 's/^.*;MQ=\([0-9]*\.[0-9]*\);.*$/\1/' > ${file}_MQ.txt

done
