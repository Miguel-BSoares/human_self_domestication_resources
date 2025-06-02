##################################################################################
##################################################################################
## The following scripts follow GATK best practices for germline variant discovery
## I deposited a simplified version that shows the skeleton and recipe of the pipeline
## Please adapt it to your needs and feel free to omit steps that are not necessarity for your research
## I am using <reference> to the refence alignment file and using lists of samples as example - both can be edited
## The <reference> in this case is a pseudo-reference


####################################################
## REMOVE READ SMALL READ LENGHTS BEFORE ALIGNING ##
####################################################

for file in *fastq.gz
        do seqkit seq -m 35 $file -g |
        gzip > seqkit_filter/${file}_35.fastq.gz

##
for sample in $list
        do      bwa aln -n 0.03 -l 1024 <reference> ${sample}.fastq.gz |
                bwa samse <reference> - ${sample}.sam ${sample}.fastq.gz |
                samtools view -b -F 4 |
                samtools sort --threads 10 > bwa_aln/${sample}.bam
done

#####################
## make dictionary ##
#####################

picard32 CreateSequenceDictionary \ 
      R=<reference> \ 
      O=<reference>.dict

######################
## revert sam files ##
######################

list="SAMPLE_A
SAMPLE_B"

for sample in $list

do picard32 RevertSam \
        I=${sample}.bam \
        O=unmapped_bam/${sample}_revert_sam.bam \
        SANITIZE=TRUE \
        ATTRIBUTE_TO_CLEAR=XT \
        ATTRIBUTE_TO_CLEAR=XN \
        ATTRIBUTE_TO_CLEAR=AS \
        ATTRIBUTE_TO_CLEAR=OC \
        ATTRIBUTE_TO_CLEAR=OP \
        TMP_DIR=tmp/ \

done

######################
  ## merge bams ##
######################

list="SAMPLE_A
SAMPLE_B"

for sample in $list

do picard MergeBamAlignment --ALIGNED_BAM ${sample}.bam \
        -UNMAPPED_BAM unmapped_bam/${sample}_revert_sam.bam -O merged_alignments/${sample}.merged.bam \
        -VALIDATION_STRINGENCY SILENT \
        -R <reference> \
        --TMP_DIR tmp/ \

done

################################
## mark and remove duplicates ##
################################

list="SAMPLE_A
SAMPLE_B"

for sample in $list

do picard MarkDuplicates \
        --INPUT replaced_readgroup/${sample}.bam \
        --OUTPUT removed_duplicates/${sample}.sorted.removed.bam \
        --METRICS_FILE dup_metrics.txt \
        --MAX_RECORDS_IN_RAM 99999 \
        --REMOVE_DUPLICATES true \
        --TMP_DIR tmp/ \

done

################################
    ## haplotype calling ##
################################

for sample in $list

        do gatk HaplotypeCaller \
        -R <reference> \
        -I removed_duplicates/${sample}.sorted.removed.bam \
        -O gvfc/${sample}.g.vcf.gz \
        -ERC GVCF \
        --min-base-quality-score 30 \
        --tmp-dir tmp/ \

done

################################
 ## rename the samples in vcf ##
################################

list="SAMPLE_A
SAMPLE_B"

for sample in $list

        do picard RenameSampleInVcf -I ${sample}.g.vcf.gz -O renamed_vcf/${sample}.g.vcf.gz --NEW_SAMPLE_NAME ${sample}

done

##################################
 ## combine the individual vcfs ##
##################################
#can be looped#

gatk CombineGVCFs -R ../../random_database/random_db.fasta \
        --variant ${sample}.g.vcf.gz \
        --variant ${sample}.g.vcf.gz \
        --variant ${sample}.g.vcf.gz \
        --variant ${sample}.g.vcf.gz \
        -O combined/human_random.g.vcf.gz --tmp-dir ../tmp/

##################################
      ## Select just SNPs ##
##################################

##################################
     ## variant filtration ##
##################################

gatk VariantFiltration \
                -R <reference> \
                -V human_snps.vcf.gz \
                --window 35 \
                --cluster 3 \
				-filter "QD < 2.0" --filter-name "QD2" \
                -filter "QUAL < 30.0" --filter-name "QUAL30" \
                -filter "SOR > 2.0" --filter-name "SOR2" \
                -filter "FS > 60.0" --filter-name "FS60" \
                -filter "MQ < 60.0" --filter-name "MQ60" \
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \ 
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                -O human_snps_filtered.vcf.gz



#######################
##bcf filters post GATK
#count positions in vcf/bcf files - to be done after filtering steps

bcftools query -f '%POS\n' $filename_filtered.vcf | wc -l

#show flags for each filter 
zgrep -v "^#" snps_filtered_v5.vcf.gz | cut -f 7 | sort -T tmp/| uniq -c	

#collect only those who PASS filters
bcftools view -f 'PASS' human_snps_filtered.vcf.gz > human_snps_filtered_pass.vcf

#filter for bialellic alleles, missing data and minor allelic frequency
bcftools view --max-alleles 2 -i 'F_MISSING<0.1' -q 0.05:minor human_snps_filtered_pass.vcf >  human_unique_biallelic_g01_maf005.vcf

#filter for trialellic alleles, missing data and minor allelic frequency
bcftools view --max-alleles 3 -i 'F_MISSING<0.1' -q 0.05:minor human_snps_filtered_pass.vcf >  human_unique_triallelic_g01_maf005.vcf

#filter for different average allelic depths
bcftools view  -i  'AVG(FMT/AD)>30' human_unique_biallelic_g01_maf005.vcf > human_unique_biallelic_g01_maf005_AD30.vcf
bcftools view  -i  'AVG(FMT/AD)>10' human_unique_biallelic_g01_maf005.vcf > human_unique_biallelic_g01_maf005_AD10.vcf

bcftools view  -i  'AVG(FMT/AD)>30' human_unique_triallelic_g01_maf005.vcf > human_unique_triallelic_g01_maf005_AD30.vcf
bcftools view  -i  'AVG(FMT/AD)>10' human_unique_triallelic_g01_maf005.vcf > human_unique_triallelic_g01_maf005_AD10.vcf
