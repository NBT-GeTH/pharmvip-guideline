#!/bin/bash
# while getopts u:a:f: flag
# do
#     case "${flag}" in
#         u) username=${OPTARG};;
#         a) age=${OPTARG};;
#         f) fullname=${OPTARG};;
#     esac
# done

# FILES="./data/*.vcf"
# # FILES2="./data/tester.vcf"
# for f in $FILES
# do
#     echo "Processing ${f} file..."
#     bgzip < ${f} > "${f}.gz"
#     tabix -p vcf "${f}.gz" 
# done

# #!/bin/bash
# FILES="./bigchunk/*.vcf.gz"
# for f in $FILES
# do
#     echo "Processing ${f} file..."
#     gunzip -k -f ${f}
# done

#!/bin/bash
FILES="./bigchunk/*.vcf"
for f in $FILES
do
    echo "Processing ${f} file..."
    gzip -k  ${f}
done