# PharmVIP - Guideline Package

PharmVIP (Pharmacogenomic Variant Analysis and Interpretation Platform) Guideline Package.
This package offers allelic determination of 17 pharmacogenes from the Next-generation sequencing (NGS) genotypic data and then reports the relevant Clinical Pharmacogenetics Implementation Consortium (CPIC) drug guideline recommendations based on the predicted allele.

## Setup

### Dependencies
*   Python 3.8+.
*   [pandas 1.2.4](https://pandas.pydata.org/)
*   [Setuptools 56.2.0](https://setuptools.readthedocs.io/en/latest/)
*   [cyvcf2 0.30.4](https://github.com/brentp/cyvcf2)

You can set up Python virtual environment (you might need to install the
`python3-venv` package first) with all needed dependencies inside the repository using:

```shell
python3 -m venv pharmvip
source pharmvip/bin/activate
pip install -r requirements.txt 
```

### Install as a python package

This package can be installed as a python package by following the below instruction: 

On the repository root directory, run this command to install the pacakage.
```shell
pip install -e .
```

## Usage 

We provide a bash script for running the package, which is located at `/resources/scripts/pharmvip_guideline_allele_matcher.sh`
edit file to set the input path and run by :
```shell
bash /resources/scripts/pharmvip_guideline_allele_matcher.sh
```

### Script fomat

```shell
pharmvip_guideline allele_matcher \
    --ana_user_id 1 \
    --ana_id 1 \
    --ana_options_cpic true \
    --ana_best_candidate true \
    --ana_genes_cyp2d6 true \
    --ana_options_hla true \
    --vcf_gz_file ${PATH_TO_INPUT_FILE} \
    --diplotype_cyp2d6 "$SCRIPTPATH/../samples/HS01011/diplotype_CYP2D6.tsv" \
    --diplotype_hla "$SCRIPTPATH/../samples/HS01011/diplotype_HLA.tsv" \
    --outputs ${PATH_TO_OUTPUT_FILE} \
    --dbpmcgenomics ${PATH_TO_OUTPUT_FILE}
```
## Input requirement

To ensure that the package run properly. The input file should be processed as described below.

### Default input GVCF of PharmVIP 
The GVCF should look like this:
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
chr1	97077743	.	T	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077744	.	T	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077745	.	T	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077746	.	A	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077747	.	A	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077748	.	A	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077749	.	A	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077750	.	T	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077751	.	G	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077752	.	C	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077753	.	T	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
chr1	97077754	.	T	.	135.82	.	DP=36	GT:AD:DP:RGQ	0/0:36,0:36:99
```

### Input preprocessing

Input preprocessing can be descript as follow :

<!-- ![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true) -->
![alt text](https://github.com/NBT-GeTH/pharmvip-guideline/blob/master/resources/samples/vcf_processing.png )

#### 1. Map to Reference Grch38
```shell
bwa mem -K 10000000 \
       -M reference.fasta\
        -R '@RG\\tID:input\\tLB:input\\tSM:input\\tPL:ILLUMINA' \
        $read1 \
        $read2 | samtools view -bS -@ ${task.cpus} -o input.bam
```
#### 2. Mark Duplicates
```shell
java -jar picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt
```
#### 3. Base Recalibration
```shell
gatk BaseRecalibrator \
   -I input.bam \
   -R reference.fasta \
   --known-sites sites_of_variation.vcf \
   --known-sites another/optional/setOfSitesToMask.vcf \
   -O recal_data.table
```
```shell
gatk PrintReads \
	-R reference.fasta \
   	-I input.bam \
	-BQSR recal_data.table \
	-o input.recal.bam
```
```shell
samtools index input.recal.bam input.recal.bam.bai
```    

#### 4. Call Variants Per-Sample
```shell
gatk --java-options "-Xmx4g" HaplotypeCaller \
-R reference.fasta \
-I input.bam \
-ERC GVCF \
-variant_index_type LINEAR  \
-variant_index_parameter 128000 \
-o input.g.vcf
```

#### 5. Joint-Call Reference
```shell
gatk --java-options "-Xmx4g" GenotypeGVCFs \
-R reference.fasta \
--include-non-variant-sites true \
--allow-old-rms-mapping-quality-annotation-data true \
-V input.g.vcf.gz \
-O input.vcf.gz \
```

#### 6. Consolidate GVCFs
```shell
gatk --java-options "-Xmx4g" CombineGVCFs \
-R reference.fasta \
--variant sample1.g.vcf.gz \
--variant sample2.g.vcf.gz \
-O combineCohort.g.vcf.gz
```

#### 7. Joint-Call Cohort
```shell
gatk --java-options "-Xmx4g" GenotypeGVCFs \
-R reference.fasta \
--include-non-variant-sites true \
--allow-old-rms-mapping-quality-annotation-data true \
-V combineCohort.g.vcf.gz \
-O GenotypeCohort.vcf.gz
```

#### 8. Select 1 sample (GVCF)
```shell
gatk --java-options "-Xmx4g" SelectVariants \
-V GenotypeCohort.vcf.gz \
-sn sample1 \
-O sample1.vcf.gz
```

## :warning: DISCLAIMER:

PharmVIP platform is currently provided for research and informational purposes only and is not intended to
diagnose, treat, cure, or prevent any disease. The report generated by the platform is not meant to be a
substitute for professional medical advice, diagnosis, or treatment provided by a physician or other qualified
health care professional. Only a physician, pharmacist or other healthcare professional should advise a
patient on the use of the medications prescribed. Any application of the content provided or obtained through
the use of the platform is therefore solely at the user’s own risk and responsibility.

## Contact

   :world_map: [National Biobank of Thailand (NBT)](https://goo.gl/maps/PUMwh6WKUvGNJeym7)<br />
   :phone: Tel : 025647000 Ext. 71475<br />
   :email: nbt@nstda.or.th<br />
