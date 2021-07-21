#!/bin/bash
# while getopts u:a:f: flag
# do
#     case "${flag}" in
#         u) username=${OPTARG};;
#         a) age=${OPTARG};;
#         f) fullname=${OPTARG};;
#     esac
# done

# FILES="../samples/1000t/*.vcf"
# # FILES2="./data/tester.vcf"
# for f in $FILES
# do
#     echo "Processing ${f} file..."
#     # bgzip < ${f} > "${f}.gz"
#     # gunzip -k -f ${f}
#     # tabix -p vcf "${f}.gz" 
# done

# FILES="../samples/1000t/*.vcf.gz"
# for f in $FILES
# do
#     SUBSTRING=$(echo  ${f} | rev | cut -c4- | rev)
#     echo $SUBSTRING
#     if test -f $SUBSTRING; then
#         echo "${SUBSTRING} exists."
#     else 
#         gunzip -k -f $SUBSTRING
#         echo "${SUBSTRING} does not exist."
#     fi
# done


FILES="../samples/1000t/*.vcf"
for f in $FILES
do
    SUBSTRING=$(echo  "${f}.gz.tbi")
    echo $SUBSTRING
    if test -f $SUBSTRING; then
        echo "${SUBSTRING} exists."
    else 
        echo "${SUBSTRING} does not exist."
        bgzip < ${f} > "${f}.gz"
        tabix -p vcf "${f}.gz" 

    fi
    # bgzip < ${f} > "${f}.gz"
    # gunzip -k -f ${f}
    # tabix -p vcf "${f}.gz" 
done

#!/bin/bash
# FILES="./bigchunk/*.vcf"
# for f in $FILES
# do
#     echo "Processing ${f} file..."
#     gzip -k  ${f}
# done