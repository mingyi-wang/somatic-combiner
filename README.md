# somatic-combiner

Usage:
java -jar somaticCombiner.jar -v Vardict.vcf -l Lofreq_snvs.vcf -L Lofreq_indels.vcf -u Muse.vcf -M Mutect2.vcf.gz -m Mutect.vcf -s strelka.snvs.vcf -S strelka.indels.vcf -o output.vcf

E.g.
java -jar somaticCombiner.jar -v T:\DCEG\CGF\Bioinformatics\Production\Mingyi\Vardict\dream_set4_fixed_vt_sorted_psssed.vcf -l T:\DCEG\Projects\Exome\builds\build_Dream\Results\Lofreq\Lofreq_set4_WGS_somatic_final_minus-dbsnp.snvs.vcf.gz -L T:\DCEG\Projects\Exome\builds\build_Dream\Results\Lofreq\Lofreq_set4_WGS_somatic_final_minus-dbsnp.indels_vt_sorted.vcf.gz -u T:\DCEG\Projects\Exome\builds\build_Dream\Results\Muse\set4.vcf.gz -M T:\DCEG\Projects\Exome\builds\build_Dream\Results\Mutect2\merged_set4_all_raw_fixed_sorted_vt_sorted.vcf.gz -m T:\DCEG\Projects\Exome\builds\build_Dream\Results\Mutect\merged_set4_all.vcf -s T:\DCEG\Projects\Exome\builds\build_Dream\Results\Strelka\set4\results\variants\somatic.snvs.vcf -S T:\DCEG\Projects\Exome\builds\build_Dream\Results\Strelka\set4\results\variants\somatic.indels_vt_sorted.vcf.gz -o T:\DCEG\Home\wangm6\tmp2\merged.vcf

1. All the VCFs with INDELs must vt normalization
2. The program can merge any two or more VCFs.	
3. For the individual VCF, in the header line, the sample columns must contain "TUMOR" or "NORMAL" (case insenstitive though) to distiguish two samples. The Lofreq VCF can leave empty since it doesn't provide sample columns.
3. "GT:DP:AD" values in the output VCF will take the values according to priority order: Mutect2 > Mutect > Muse > Vardict > Strelka > Lofreq.
4. The output VCF is a superset of all VCFs. The variants calling status is tagged by "NumCallers=6;lLsSumMv=10101111" in INFO column. 