#!/bin/bash
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Add description of the script functions here."
   echo
   echo "Syntax: scriptTemplate [-g|d] [PATHtoDir|PATHtoFile]"
   echo "options:"
   echo "-h,--help      Print Help."
   echo "-d             decompose .vcf.gz."
   echo "-i|--indexe    from vcf create vcf.gz and vcf.gz.tbi"
   echo
}

Ungz() {
    extension="${INPP##*.}"
    if (( $extension  == "gz" ))
    then 
        FILES="${INPP}"
    else
        FILES="${INPP}/*.vcf.gz"
    fi
    for f in $FILES
    do
        echo "Processing ${f} file..."
        SUBSTRING=$(echo  ${f} | rev | cut -c4- | rev)
        gzip  -d <  ${f} > ${SUBSTRING}
    done
}

IndexManipulate() {
    extension="${INPP##*.}"
    if (( $extension  == "vcf" ))
    then 
        FILES="${INPP}"
    else
        FILES="${INPP}/*.vcf"
    fi
    
    for f in $FILES
    do
        echo "Processing ${f} file..."
        bgzip -f < ${f} > "${f}.gz"
        tabix -p vcf "${f}.gz" 
    done
}

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    -d|--decompose)
      INPP="$2"
      shift # past argument
      shift # past value
      Ungz
      ;;
    -h|--help)
      SEARCHPATH="$2"
      shift # past argument
      shift # past value
      Help
      ;;
    -i|--indexe)
      INPP="$2"
      shift # past argument
      shift # past value
      IndexManipulate
      ;;
    --default)
      DEFAULT=YES
      shift # past argument
      ;;
    *)    # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      Help
      ;;
  esac
done
