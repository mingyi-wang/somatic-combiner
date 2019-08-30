# somatic-combiner
A consensus ensemble approach which can combine somatic variants from six popular callers: Vardict, Lofreq, Muse, Mutect2, Mutect, Strelka.

Usage:
java -jar somaticCombiner.jar -v Vardict.vcf -l Lofreq_snvs.vcf -L Lofreq_indels.vcf -u Muse.vcf -M Mutect2.vcf.gz -m Mutect.vcf -s strelka.snvs.vcf -S strelka.indels.vcf -o output.vcf

E.g.
Using example folder files to run test
java -jar somaticCombiner.jar -v example/dream_set4_fixed_vt_sorted_psssed.vcf -l example/Lofreq_set4_WGS_somatic_final_minus-dbsnp.snvs.vcf.gz -L example/Lofreq_set4_WGS_somatic_final_minus-dbsnp.indels_vt_sorted.vcf.gz -u example/set4.vcf.gz -M example/merged_set4_all_raw_fixed_sorted_vt_sorted.vcf.gz -m /example/merged_set4_all.vcf -s example/somatic.snvs.vcf -S example/somatic.indels_vt_sorted.vcf.gz -o example/merged.vcf

1. All the VCFs must called from same BAMs and reference genomes. The VCFs with INDELs must do vt normalization before merge.
2. The program can merge any two or more VCFs.
3. Only "PASS" in the FILTER column will be used for merge.
4. For the individual VCF, in the header line, the sample columns must contain "TUMOR" or "NORMAL" (case insenstitive though) to distiguish two samples. The Lofreq VCF can leave empty since it doesn't provide sample columns.
5. "GT:DP:AD" values in the output VCF will take the values according to priority order: Mutect2 > Mutect > Muse > Vardict > Strelka > Lofreq.
6. The output VCF is a superset of all VCFs. The variants calling status is tagged by "NumCallers=6;lLsSumMv=10101111" in INFO column. 