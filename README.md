# somatic-combiner v1.0
A consensus ensemble approach which can combine somatic variants VCFs generated from seven popular callers: LoFreq, MuSE, MuTect2, MuTect, Strelka, VarDict and VarScan.

## Prerequisite
java jdk or jre v1.8 or above

## Installation
git clone https://github.com/mingyi-wang/somatic-combiner.git

## Usage
- Command line:

The input VCF files can be .gz files.

java -jar somaticCombiner.jar -L ${LoFreq_INDEL_VCF} -l ${LoFreq_SNV_VCF} -u ${MuSE_VCF} -M ${MuTect2_VCF} -m ${MuTect_VCF} -s ${Strelka_SNV_VCF} -S ${Strelka_INDEL_VCF} -v ${VarScan_SNV_VCF} -V ${Varscan_INDEL_VCF} -D ${VarDict_VCF} -o ${OUTPUT_VCF}

- Example

Using the example VCF files in the example folder to run test

cd example && sh ./run_combine_example.sh

## Notes
1. All input VCFs must be called from same tumor/normal paired BAMs and reference genome. The VCFs with INDELs must do vt normalization before merge.
2. The program can merge any two or more VCFs from above seven callers.
3. Only the variants with a "PASS" in the FILTER column will be used for the merge process.
4. For the individual VCF, in the header line, the sample columns must contain "TUMOR" or "NORMAL" (case insenstitive though) to distiguish two samples. The Lofreq VCF can leave empty for those columns since it does not return sample columns.
5. "GT:DP:AD" values in the output VCF will use the values from input VCFs according to the configured priority order: MuTect2 > MuTect > MuSE > VarDict > Strelka > LoFreq > VarScan.
6. The INFO column of an output VCF contains two fields, Tumor_AD (Tumor allelic fraction) and Tumor_DP (Tumor depth) and both of them are retrieved from individual VCFs accoring to above priority order.
6. The output VCF is a superset of all VCFs. In the output VCF, the calling status will be annotated by in the INFO column (e.g., "NumCallers=7;lLsSumMDvV=1010111110") and tagged as "PASS", "LowConf", "ADJ_LowConf" or "ADJ_PASS" in the FILTER column. The description of these two items are presented in the header part of output VCFs.
