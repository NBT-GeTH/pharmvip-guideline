
#!/bin/bash
FILES="./data/*.vcf"
# FILES2="./data/tester.vcf"
for f in $FILES
do
    echo "Processing ${f} file..."
    bgzip < ${f} > "${f}.gz"
    tabix -p vcf "${f}.gz" 
done
