
#!/bin/bash
# NOTE : Quote it else use array to avoid problems #
FILES="$1/*.vcf"
# echo $FILES
# FILES2="./data/tester.vcf"
for f in $FILES
do
    echo "Processing ${f} file..."
    bgzip < ${f} > "${f}.gz"
    tabix -p vcf "${f}.gz" 
done
