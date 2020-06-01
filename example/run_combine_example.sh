#!/bin/sh
SCRIPT=$(readlink -f "$0")
SCRIPT_DIR=$(dirname "$SCRIPT")

CMD="java -jar ${SCRIPT_DIR}/somaticCombiner.jar -L ${SCRIPT_DIR}/Lofreq_indels_vt.vcf.gz -l ${SCRIPT_DIR}/Lofreq_snvs.vcf.gz -u ${SCRIPT_DIR}/Muse.vcf.gz -M ${SCRIPT_DIR}/Mutect2_vt.vcf.gz -m ${SCRIPT_DIR}/Mutect.vcf.gz -s ${SCRIPT_DIR}/Strelka_snvs.vcf.gz -S ${SCRIPT_DIR}/Strelka_indels_vt.vcf.gz -v ${SCRIPT_DIR}/Varscan_snvs.vcf.gz -V ${SCRIPT_DIR}/Varscan_indels_vt.vcf.gz -D ${SCRIPT_DIR}/Vardict_vt.vcf.gz -t 0.1 -o ${SCRIPT_DIR}/combined_7callers.vcf"
echo $CMD
eval $CMD
if [[ $? -ne 0 ]]; then
   echo "Error: combine is failed!"
   exit 1
else
   echo "Done!"
fi
