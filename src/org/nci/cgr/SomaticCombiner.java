package org.nci.cgr;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.*;
public class SomaticCombiner {
	public static ArrayList<Caller> callerList=null;
	public final static String COUNT_TAG="NumCallers";
	public static String callerSymbols="";
	public static String callerNames="";
	public static ArrayList<String> medianAFs=new ArrayList<String>();
	public static void main(String[] args) throws ClassNotFoundException, SQLException, IOException {
		// TODO Auto-generated method stub
		// TODO Auto-generated method stub
		// Parameter examples:
		// --lofreq-snv
		// --lofreq-indel T:\DCEG\Projects\Exome\builds\2014-12-17\Ensemble_New_Annotation\variants\variants_annotated_new.vcf T:\DCEG\Home\wangm6\Bin\variant_annotation_scripts\tmp.txt
		// --vardict
		// --strelka-snv
		// --strelka-indel
		// --muse
		// --mutect
		// --mutect2
		// -o 
		//
		
		// augument list:
		// -v T:\DCEG\CGF\Bioinformatics\Production\Mingyi\Vardict\dream_set4_fixed_vt_sorted_psssed.vcf 
		// -l T:\DCEG\Projects\Exome\builds\build_Dream\Results\Lofreq\Lofreq_set4_WGS_somatic_final_minus-dbsnp.snvs.vcf.gz 
		// -L T:\DCEG\Projects\Exome\builds\build_Dream\Results\Lofreq\Lofreq_set4_WGS_somatic_final_minus-dbsnp.indels_vt_sorted.vcf.gz 
		// -u T:\DCEG\Projects\Exome\builds\build_Dream\Results\Muse\set4.vcf.gz 
		// -M T:\DCEG\Projects\Exome\builds\build_Dream\Results\Mutect2\merged_set4_all_raw_fixed_sorted_vt_sorted.vcf.gz 
		// -m T:\DCEG\Projects\Exome\builds\build_Dream\Results\Mutect\merged_set4_all.vcf 
		// -s T:\DCEG\Projects\Exome\builds\build_Dream\Results\Strelka\set4\results\variants\somatic.snvs.vcf 
		// -S T:\DCEG\Projects\Exome\builds\build_Dream\Results\Strelka\set4\results\variants\somatic.indels_vt_sorted.vcf.gz -o T:\DCEG\Home\wangm6\tmp2\merged.vcf
		
		// Parameter:
		// -l T:\\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter2_downsample\Lofreq\Lofreq_NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_WES_somatic_final.snvs.vcf.gz 
		// -L T:\\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter2_downsample\Lofreq\Lofreq_NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_WES_somatic_final.indels.vcf.gz 
		// -u T:\\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter2_downsample\Muse\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_passed.vcf.gz 
		// -M T:\\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter2_downsample\Sentieon\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_WES_vt_sorted.vcf.gz 
		// -s T:\\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter2_downsample\Strelka\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40\results\variants\somatic.snvs.vcf.gz 
		// -S T:\\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter2_downsample\Strelka\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40\results\variants\somatic.indels_vt_sorted.vcf.gz 
		// -o T:\\DCEG\Home\wangm6\tmp2\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_4callers_voting.vcf
		
		
		// -l T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_production_downsample\Lofreq\Lofreq_NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_WES_somatic_final.snvs.vcf.gz -L T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_production_downsample\Lofreq\Lofreq_NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_WES_somatic_final.indels.vcf.gz -u T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_production_downsample\Muse\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_passed.vcf.gz -M T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_production_downsample\Sentieon\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_WES_vt_sorted_fixed.vcf -s T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_production_downsample\Strelka\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40\results\variants\somatic.snvs.vcf.gz -S T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_production_downsample\Strelka\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40\results\variants\somatic.indels_vt_sorted.vcf.gz -o T:\DCEG\Home\wangm6\tmp2\NA24385_tumor10PNA12878-1_1000_vs_NA24385_germline1_40_4callers_voting.vcf
		// -l T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_normal_pipeline\Lofreq\Lofreq_CTRL_NA24385_tumor50PNA12878-1_vs_CTRL_NA24385_germline1_WES_somatic_final.snvs.vcf.gz -L T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_normal_pipeline\Lofreq\Lofreq_CTRL_NA24385_tumor50PNA12878-1_vs_CTRL_NA24385_germline1_WES_somatic_final.indels.vcf.gz -u T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_normal_pipeline\Muse\CTRL_NA24385_tumor50PNA12878-1_vs_CTRL_NA24385_germline1_passed.vcf.gz -M T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_normal_pipeline\Sentieon\CTRL_NA24385_tumor50PNA12878-1_vs_CTRL_NA24385_germline1_WES_vt_sorted_fixed.vcf -s T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_normal_pipeline\Strelka\CTRL_NA24385_tumor50PNA12878-1_vs_CTRL_NA24385_germline1\results\variants\somatic.snvs.vcf.gz -S T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_normal_pipeline\Strelka\CTRL_NA24385_tumor50PNA12878-1_vs_CTRL_NA24385_germline1\results\variants\somatic.indels_vt_sorted.vcf.gz -o T:\DCEG\Home\wangm6\tmp2\CTRL_NA24385_tumor50PNA12878-1_vs_CTRL_NA24385_germline1_4callers_voting.vcf
		// -l T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Lofreq\Lofreq_NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_WES_somatic_final.snvs.vcf.gz -L T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Lofreq\Lofreq_NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_WES_somatic_final.indels.vcf.gz -u T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Muse\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_passed.vcf.gz -M T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Sentieon\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_WES_vt_sorted_fixed.vcf -s T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Strelka\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40\results\variants\somatic.snvs.vcf.gz -S T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Strelka\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40\results\variants\somatic.indels_vt_sorted.vcf.gz -o T:\DCEG\Home\wangm6\tmp2\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_4callers_voting.vcf
		
		// -l T:\\Projects\Exome\builds\build_Dream\Results\Lofreq\Lofreq_set1_WGS_somatic_final_minus-dbsnp.snvs.vcf.gz -L T:\\Projects\Exome\builds\build_Dream\Results\Lofreq\Lofreq_set1_WGS_somatic_final_minus-dbsnp.indels_vt_sorted.vcf -M T:\\Projects\Exome\builds\build_Dream\Results\Mutect2\merged_set1_all_raw_fixed_sorted_vt_sorted.vcf -s T:\\Projects\Exome\builds\build_Dream\Results\Strelka\set1\results\variants\somatic.snvs.vcf.gz -S T:\\Projects\Exome\builds\build_Dream\Results\Strelka\set1\results\variants\somatic.indels_vt_sorted.vcf.gz -u T:\\Projects\Exome\builds\build_Dream\Results\Muse\set1.vcf -o T:\\Projects\Exome\builds\build_Dream\Results\Ensemble\set1_4callers_voting_ourown.vcf
		
		// -l T:\DCEG\Projects\Exome\builds\build_Dream\Results\Lofreq\Lofreq_COLO_WGS_somatic_final_minus-dbsnp.snvs.vcf -L T:\DCEG\Projects\Exome\builds\build_Dream\Results\Lofreq\Lofreq_COLO_WGS_somatic_final_minus-dbsnp.indels_vt_sorted.vcf -M T:\DCEG\Projects\Exome\builds\build_Dream\Results\Mutect2\merged_COLO_all_raw_fixed_sorted_vt_sorted.vcf -s T:\DCEG\Projects\Exome\builds\build_Dream\Results\Strelka\COLO\results\variants\somatic.snvs.vcf.gz -S T:\DCEG\Projects\Exome\builds\build_Dream\Results\Strelka\COLO\results\variants\somatic.indels_vt_sorted.vcf.gz -u T:\DCEG\Projects\Exome\builds\build_Dream\Results\Muse\COLO.vcf -o T:\DCEG\Projects\Exome\builds\build_Dream\Results\Ensemble\COLO_4callers_voting_ourown.vcf
		
		// -l T:\DCEG\Projects\Exome\builds\build_precisionFDA_test\Results_normal_pipeline\Lofreq\Lofreq_2_80_20_WES_somatic_final.snvs.vcf.gz -L T:\DCEG\Projects\Exome\builds\build_precisionFDA_test\Results_normal_pipeline\Lofreq\Lofreq_2_80_20_WES_somatic_final.indels.vcf.gz -M T:\DCEG\Projects\Exome\builds\build_precisionFDA_test\Results_normal_pipeline\Sentieon\2_80_20_WES_vt_sorted.vcf.gz -s T:\DCEG\Projects\Exome\builds\build_precisionFDA_test\Results_normal_pipeline\Strelka\2_80_20\results\variants\somatic.snvs.vcf.gz -S T:\DCEG\Projects\Exome\builds\build_precisionFDA_test\Results_normal_pipeline\Strelka\2_80_20\results\variants\somatic.indels_vt_sorted.vcf.gz -v T:\DCEG\Projects\Exome\builds\build_precisionFDA_test\Results_normal_pipeline\Vardict_0919\Vardict_merged_2_80_20_final_vt_sorted.vcf.gz -o T:\DCEG\Projects\Exome\builds\build_precisionFDA_test\Results_normal_pipeline\Ensemble\2_80_20_4callers_voting.vcf
		initCallers();
		CommandLineParser parser=new DefaultParser();
		Options options=new Options();
		for (Caller caller :callerList)
			options.addOption(caller.getSymbol(),caller.getName(),true,caller.getDescription());
		
		options.addOption("o","output",true,"Output VCF file");
		String outputFilePath=null;
		

		CommandLine line = null;
		try{
			line=parser.parse(options, args);
		}
		catch (ParseException e) {
            System.out.println(e.getMessage());
            //formatter.printHelp("utility-name", options);

            System.exit(1);
        }
		int count=0;
		for (Caller caller:callerList) {
			if (line.hasOption(caller.getSymbol())) {
				caller.setFilePath(line.getOptionValue(caller.getSymbol()));
				File tmpFile=new File(caller.getFilePath());
				if (tmpFile.exists()) 
					count++;
				else {
					System.out.println("Error:"+caller.getFilePath()+" does not existed! Please check for that! The program will skip this VCF file!");
					System.exit(1);
				}				
			}
		}
		if (count<2) {
			System.out.println("Error: Only one VCF is available. No merge needed!");
			System.exit(1);
		}
		outputFilePath=line.getOptionValue("o");
		
		int snvCallerNum=0;
		int indelCallerNum=0;
		List<Variant> list =new ArrayList<Variant>();
		System.out.println("Starting loading ...");
		for (Caller caller:callerList) {
			if (caller.getFilePath()!=null) {
				VCFFile vcfFile=new VCFFile(caller.getFilePath(),caller.getName(),caller.getPriority(),caller.getSet());
//				if (caller.getName().contains("strelka-indel")) {
//					System.out.println("Found!");
//				}
			    if(vcfFile.importVariants(list)>0) {
			    	caller.setVcfFile(vcfFile);
			    	callerSymbols+=caller.getSymbol();
			    	callerNames+=caller.getName()+",";
			    	if (caller.getType().equals("BOTH")){
			    		snvCallerNum++;
			    		indelCallerNum++;
			    		
			    	}
			    	else 
			    		if (caller.getType().equals("SNV"))
			    			snvCallerNum++;
				    	else
				    		indelCallerNum++;
			    }
			    
			}
		}
		if (callerSymbols.length()<=1) {
			System.out.println("Caller number is less than 2. Exit!");
            System.exit(1);
		}
		callerNames=callerNames.substring(0, callerNames.length()-1);
		List<MergedVariant> mergedList=new ArrayList<MergedVariant>();
		System.out.println("Starting merging...");
		long i=0;
		for (Variant p : list) {
			i++;
			System.out.println(i);
//			if (p.getVariantContext().getStart()==6263642)
//				System.out.println(p.getVariantContext().getContig()+"\t"+p.getVariantContext().getStart());
			MergedVariant mp=new MergedVariant(p.getVariantContext(), p.getCaller(),p.getPriority(),p.getSet());
			int index=mergedList.indexOf(mp);
			if (index!=-1) {
				MergedVariant mergedVariant=mergedList.get(index).merge(mp);
				if (mergedVariant==null)
					mergedList.add(new MergedVariant(p.getVariantContext(), p.getCaller(),p.getPriority(),p.getSet()));
				else
				    mergedList.set(index, mergedVariant);
			}
			else {
				
				mergedList.add(mp);
			}
		}
		System.out.println("Sorting merged variants ...");
		Collections.sort(mergedList,Variant.VariantComparator);
		System.out.println("Writing VCF ...");
		
//		File fw=new File(outputFilePath);
//		OutputStream bw=new FileOutputStream(fw);
//		VCFWriter vcfWriter=new VCFWriter(fw, bw, null,true,false,true,false);
//        i=0;
//   		for (Variant p : mergedList) {
//    			i++;
//    			System.out.println("Writing:"+i);
//    			vcfWriter.add(p.getVariantContext());
//   		}     
//		vcfWriter.close();
		
		BufferedWriter bw=new BufferedWriter(new FileWriter(new File(outputFilePath)));
		VCFFile firstVCFFile=prepareHeader();	
		VCFFile.writeHeader(firstVCFFile.getmHeader(), bw);
		
//		System.out.println("Calculating Median AF for each fragment ...");
//		String currentContig="";
//		int currentStart=0;
//		i=0;
//		ArrayList<Float> tumorAFs=new ArrayList<Float>();
//		for (Variant p : mergedList) {
//			if (i==0) {
//				currentContig=p.variantContext.getContig();
//			    currentStart=p.variantContext.getStart();
//			}
//
//
//			if ((p.variantContext.getContig()==currentContig)&&(p.variantContext.getStart()<=currentStart+500)){
//			
//
//			   float tumorAF=firstVCFFile.segmentVariants(p,snvCallerNum,indelCallerNum);
//			   if (tumorAF>0 && tumorAF<=1)
//				   tumorAFs.add(tumorAF);
//			}
//			else {
//				if (tumorAFs.size()>0) {
//					Collections.sort(tumorAFs);
//					float median;
//					if(tumorAFs.size()%2==0) {
//				      float sumOfMiddleElements= tumorAFs.get(tumorAFs.size()/2)+tumorAFs.get(tumorAFs.size()/2-1);
//				      median=((float)sumOfMiddleElements/2);
//					}
//					else
//						median=(float)tumorAFs.get(tumorAFs.size()/2);
//					if (median<=0.1) {
//						medianAFs.add(currentContig+":"+currentStart+":"+median+":"+tumorAFs.size());
//						System.out.println(currentContig+":"+currentStart+":"+median+":"+tumorAFs.size());
//					}
//					tumorAFs.removeAll(tumorAFs);
//				}
//				if (!p.variantContext.getContig().equals(currentContig)) {
//					currentContig=p.variantContext.getContig();
//				    currentStart=p.variantContext.getStart();
//				}
//				else {
//					currentStart=p.variantContext.getStart();					
//				}
//				float tumorAF=firstVCFFile.segmentVariants(p,snvCallerNum,indelCallerNum);
//				if (tumorAF>0 && tumorAF<=1)
//					   tumorAFs.add(tumorAF);
//			}
//			i++;
//		}
		
		
		i=0;
		for (Variant p : mergedList) {
			i++;
			System.out.println("Writing:"+i);
//			if (i==3439)
//				System.out.println("found!");
			firstVCFFile.writeVariants(bw,p,snvCallerNum,indelCallerNum);
//			firstVCFFile.writeVariants(bw,p,snvCallerNum,indelCallerNum,medianAFs);
		}
		bw.close();
		System.out.println("Done!");
	}
	
	private  static  void initCallers() {
		SomaticCombiner.callerList=new ArrayList<Caller>();
		callerList.add(new Caller("lofreq-snv", "Lofreq SNV VCF file" , "l", (byte) 0b0010000, 2,"SNV","Lofreq"));
		callerList.add(new Caller("lofreq-indel", "Lofreq INDEL VCF file" , "L", (byte) 0b0010000, 2,"INDEL","Lofreq"));
		callerList.add(new Caller("strelka-snv", "Strelka SNV VCF file" , "s", (byte) 0b0001000, 3,"SNV","Strelka"));
		callerList.add(new Caller("strelka-indel", "Strelka INDEL VCF file" , "S", (byte) 0b0001000, 3,"INDEL","Strelka"));
		callerList.add(new Caller("muse", "Muse VCF file" , "u", (byte) 0b0000100, 5,"SNV","Muse"));
		callerList.add(new Caller("mutect", "Mutect VCF file" , "m", (byte) 0b1000000, 6,"SNV","Mutect"));
		callerList.add(new Caller("mutect2", "Mutect2 VCF file" , "M", (byte) 0b0100000, 7,"BOTH","Mutect2"));
		callerList.add(new Caller("vardict", "Vardict VCF file" , "D", (byte) 0b0000010, 4,"BOTH","Vardict"));
		callerList.add(new Caller("varscan-snv", "Varscan SNV VCF file" , "v", (byte) 0b0000001, 1,"SNV","Varscan"));
		callerList.add(new Caller("varscan-indel", "Varscan INDEL VCF file" , "V", (byte) 0b0000001, 1,"INDEL","Varscan"));
	}

	public static String callerName(String name) {		
		for (Caller caller:callerList) {
			if (name.equals(caller.getName()))
			return caller.getCallerName();
			
		}
		return "";
	}
	
	public static String nameFromSymbol(String s) {		
		for (Caller caller:callerList) {
			if (s.equals(caller.getSymbol()))
			return caller.getName();
			
		}
		return "";
	}
	
	
	public static String callerNameFromPriority(int p) {
		for (Caller caller:callerList) {
			if (p==caller.getPriority())
			return caller.getCallerName();
			
		}
		return "";	
	}
	private static VCFFile prepareHeader() {
		
		boolean firstVCF=true;
		VCFFile firstVCFFile=null;
		for (Caller caller:callerList) {
			if (caller.getVcfFile()!=null) {
				if (firstVCF) {
					firstVCFFile=caller.getVcfFile();					
					Set<VCFHeaderLine> mInfoMetaData = new HashSet<VCFHeaderLine>();
					VCFHeaderLine vcfInfoHeaderLine=new VCFInfoHeaderLine(COUNT_TAG, 1, VCFHeaderLineType.Integer,"The number of callers" );
					mInfoMetaData.add(vcfInfoHeaderLine);
					vcfInfoHeaderLine=new VCFInfoHeaderLine(callerSymbols,1, VCFHeaderLineType.String,"Calling descision of the callers in a binary string: "+callerNames );
					
					mInfoMetaData.add(vcfInfoHeaderLine);
					if (callerNames.contains("lofreq")) {
						vcfInfoHeaderLine=new VCFInfoHeaderLine("Lofreq_QUAL", 1, VCFHeaderLineType.Float,"Lofreq QUAL" );
						mInfoMetaData.add(vcfInfoHeaderLine);
					}
					if (callerNames.contains("vardict")) {
						vcfInfoHeaderLine=new VCFInfoHeaderLine("Vardict_QUAL", 1, VCFHeaderLineType.Float,"Vardict QUAL" );						
						mInfoMetaData.add(vcfInfoHeaderLine);
					}
					vcfInfoHeaderLine=new VCFInfoHeaderLine("Tumor_AF", 1, VCFHeaderLineType.Float,"Tumor Allelic Fraction from individual VCFs based on callers priority" );						
					mInfoMetaData.add(vcfInfoHeaderLine);
					vcfInfoHeaderLine=new VCFInfoHeaderLine("Tumor_DP", 1, VCFHeaderLineType.Integer,"Tumor Depth retrieved from individual VCFs based on callers priority" );						
					mInfoMetaData.add(vcfInfoHeaderLine);
					VCFHeaderLine vcfFormatHeaderLine=new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String,"Genotype" );
					mInfoMetaData.add(vcfFormatHeaderLine);
					vcfFormatHeaderLine=new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer,"Total Depth" );
					mInfoMetaData.add(vcfFormatHeaderLine);
					VCFFilterHeaderLine vcfFilterHeaderLine=new VCFFilterHeaderLine("LowConf","Low confidence call");
					mInfoMetaData.add(vcfFilterHeaderLine);
					vcfFilterHeaderLine=new VCFFilterHeaderLine("PASS","high confidence call");
					mInfoMetaData.add(vcfFilterHeaderLine);
					vcfFilterHeaderLine=new VCFFilterHeaderLine("ADJ_PASS","Adjusted high confidence call");
					mInfoMetaData.add(vcfFilterHeaderLine);
					vcfFilterHeaderLine=new VCFFilterHeaderLine("ADJ_LowConf","Adjusted low confidence call");
					mInfoMetaData.add(vcfFilterHeaderLine);
					firstVCFFile.mergeHeader(new VCFHeader(mInfoMetaData));					
					firstVCF=false;
				}
				else
					firstVCFFile.mergeHeader(caller.getVcfFile().getmHeader());
			}
		}		
		return firstVCFFile;
		
	
	}
}
