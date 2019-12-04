#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")
# zcat Lofreq_set4_WGS_somatic_final_minus-dbsnp.indels_vt_sorted.vcf.gz | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Lofreq_indels.vcf
# zcat Lofreq_set4_WGS_somatic_final_minus-dbsnp.snvs.vcf.gz | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Lofreq_snvs.vcf
# zcat muse_set4.vcf.gz | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Muse.vcf
# zcat mutect2_set4.vcf.gz | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Mutect2.vcf
# cat mutect_set4.vcf | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Mutect.vcf
# zcat strelka_set4_snvs.vcf.gz | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Strelka_snvs.vcf
# zcat strelka_set4_indels.vcf.gz | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Strelka_indels.vcf
# cat varscan_snps.vcf | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Varscan_snvs.vcf
# zcat varscan_indels.vcf.gz | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Varscan_indels.vcf
# cat dream_set4_fixed_vt_sorted_psssed.vcf | awk -F"\t" '$1~/^#/ || $1=="1" {print}' > Vardict.vcf

# java -jar ~/dev/somaticseq2/somaticCombiner.jar -L Lofreq_set4_WGS_somatic_final_minus-dbsnp.indels_vt_sorted.vcf.gz -l Lofreq_set4_WGS_somatic_final_minus-dbsnp.snvs.vcf.gz -u muse_set4.vcf.gz -M mutect2_set4.vcf.gz -m mutect_set4.vcf -s strelka_set4_snvs.vcf.gz -S strelka_set4_indels.vcf.gz -v varscan_snps.vcf -V varscan_indels.vcf.gz -D Vardict_merged_set4.vcf.gz -o 7callers.vcf
REFERENCE_GENOME=/DCEG/Projects/PopulationExome/EAGLE/b37/human_g1k_v37_no_decoy.fasta
module load vt jdk/1.8.0_111 bgzip tabix
vt_normalize()
{
   VCF=$1
   VCF_DIR=`dirname $VCF`
   if [[ $VCF != *".gz" ]]; then
     BASE_NAME=`basename $VCF .vcf`
   else
     BASE_NAME=`basename $VCF .vcf.gz`
     zcat $VCF > ${VCF_DIR}/${BASE_NAME}.vcf
     VCF=${VCF_DIR}/${BASE_NAME}.vcf
   fi
   REF=$2
   CMD="vt normalize -r $REFERENCE_GENOME -o ${VCF_DIR}/${BASE_NAME}_vt.vcf $VCF -n"
   echo $CMD
   eval $CMD
   if [[ $? -ne 0 ]]; then
     echo "Error: vt is failed!"
     exit 1
   fi

}

zip_vcf()
{
  VCF=$1
  VCF=$1
  VCF_DIR=`dirname $VCF`
  BASE_NAME=`basename $VCF .vcf`
   CMD="bgzip -c ${VCF_DIR}/${BASE_NAME}.vcf > ${VCF_DIR}/${BASE_NAME}.vcf.gz"
   echo $CMD
   eval $CMD
   if [[ $? -ne 0 ]]; then
     echo "Error: vt is failed!"
     exit 1
   fi
   CMD="tabix -p vcf ${VCF_DIR}/${BASE_NAME}.vcf.gz"
   echo $CMD
   eval $CMD

 
}
vt_normalize ${SCRIPT_DIR}/Lofreq_indels.vcf $REFERENCE_GENOME
zip_vcf ${SCRIPT_DIR}/Lofreq_indels_vt.vcf
vt_normalize ${SCRIPT_DIR}/Mutect2.vcf $REFERENCE_GENOME
zip_vcf ${SCRIPT_DIR}/Mutect2_vt.vcf
vt_normalize ${SCRIPT_DIR}/Strelka_indels.vcf $REFERENCE_GENOME
zip_vcf ${SCRIPT_DIR}/Strelka_indels_vt.vcf
vt_normalize ${SCRIPT_DIR}/Varscan_indels.vcf $REFERENCE_GENOME

zip_vcf ${SCRIPT_DIR}/Varscan_indels_vt.vcf
vt_normalize ${SCRIPT_DIR}/Vardict.vcf $REFERENCE_GENOME
zip_vcf ${SCRIPT_DIR}/Vardict_vt.vcf
zip_vcf ${SCRIPT_DIR}/Lofreq_snvs.vcf
zip_vcf ${SCRIPT_DIR}/Muse.vcf
zip_vcf ${SCRIPT_DIR}/Mutect.vcf
zip_vcf ${SCRIPT_DIR}/Strelka_snvs.vcf
zip_vcf ${SCRIPT_DIR}/Varscan_snvs.vcf


CMD="java -jar ${SCRIPT_DIR}/somaticCombiner.jar -L ${SCRIPT_DIR}/Lofreq_indels_vt.vcf.gz -l ${SCRIPT_DIR}/Lofreq_snvs.vcf.gz -u ${SCRIPT_DIR}/Muse.vcf.gz -M ${SCRIPT_DIR}/Mutect2_vt.vcf.gz -m ${SCRIPT_DIR}/Mutect.vcf.gz -s ${SCRIPT_DIR}/Strelka_snvs.vcf.gz -S ${SCRIPT_DIR}/Strelka_indels_vt.vcf.gz -v ${SCRIPT_DIR}/Varscan_snvs.vcf.gz -V ${SCRIPT_DIR}/Varscan_indels_vt.vcf.gz -D ${SCRIPT_DIR}/Vardict_vt.vcf.gz -o ${SCRIPT_DIR}/combined_7callers.vcf"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: combine is failed!"
   exit 1
else
   echo "Done!"
fi
