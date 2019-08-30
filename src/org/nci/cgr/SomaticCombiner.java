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
		
		// -l T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Lofreq\Lofreq_NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_WES_somatic_final.snvs.vcf.gz -L T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Lofreq\Lofreq_NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_WES_somatic_final.indels.vcf.gz -u T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Muse\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_passed.vcf.gz -M T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Sentieon\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_WES_vt_sorted_fixed.vcf -s T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Strelka\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40\results\variants\somatic.snvs.vcf.gz -S T:\DCEG\Projects\Exome\builds\build_UMI_NP0084_22047\Results_filter1_downsample\Strelka\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40\results\variants\somatic.indels_vt_sorted.vcf.gz -o T:\DCEG\Home\wangm6\tmp2\NA24385_tumor2PNA12878-1_100_vs_NA24385_germline1_40_4callers_voting.vcf
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
		for (Caller caller:callerList) {
			if (line.hasOption(caller.getSymbol()))
				caller.setFilePath(line.getOptionValue(caller.getSymbol()));		
		}
		outputFilePath=line.getOptionValue("o");
		
		int snvCallerNum=0;
		int indelCallerNum=0;
		List<Variant> list =new ArrayList<Variant>();
		System.out.println("Starting loading ...");
		for (Caller caller:callerList) {
			if (caller.getFilePath()!=null) {
				VCFFile vcfFile=new VCFFile(caller.getFilePath(),caller.getName(),caller.getPriority(),caller.getSet());
				if (caller.getName().contains("strelka-indel")) {
					System.out.println("Found!");
				}
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
			if (p.getVariantContext().getStart()==1770788)
				System.out.println("found!");
			MergedVariant mp=new MergedVariant(p.getVariantContext(), p.getCaller(),p.getPriority(),p.getSet());
			int index=mergedList.indexOf(mp);
			if (index!=-1) {
				MergedVariant mergedVariant=mergedList.get(index).merge(mp);
				if (mergedVariant==null)
					mergedList.add((MergedVariant) p);
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
		i=0;
		for (Variant p : mergedList) {
			i++;
			System.out.println("Writing:"+i);
			if (i==45)
				System.out.println("found!");
			firstVCFFile.writeVariants(bw,p,snvCallerNum,indelCallerNum);
		}
		bw.close();
		
	}
	
	private  static  void initCallers() {
		SomaticCombiner.callerList=new ArrayList<Caller>();
		callerList.add(new Caller("lofreq-snv", "Lofreq SNV VCF file" , "l", (byte) 0b001000, 1,"SNV","Lofreq"));
		callerList.add(new Caller("lofreq-indel", "Lofreq INDEL VCF file" , "L", (byte) 0b001000, 1,"INDEL","Lofreq"));
		callerList.add(new Caller("strelka-snv", "Strelka SNV VCF file" , "s", (byte) 0b000100, 2,"SNV","Strelka"));
		callerList.add(new Caller("strelka-indel", "Strelka INDEL VCF file" , "S", (byte) 0b000100, 2,"INDEL","Strelka"));
		callerList.add(new Caller("muse", "Muse VCF file" , "u", (byte) 0b000010, 4,"SNV","Muse"));
		callerList.add(new Caller("mutect", "Mutect VCF file" , "m", (byte) 0b100000, 5,"SNV","Mutect"));
		callerList.add(new Caller("mutect2", "Mutect2 VCF file" , "M", (byte) 0b010000, 6,"BOTH","Mutect2"));
		callerList.add(new Caller("vardict", "Vardict VCF file" , "v", (byte) 0b000001, 3,"BOTH","Vardict"));
		
	}

	public static String callerName(String name) {		
		for (Caller caller:callerList) {
			if (name.equals(caller.getName()))
			return caller.getCallerName();
			
		}
		return "";
	}
	
	public static String callerNameFromSymbol(String s) {		
		for (Caller caller:callerList) {
			if (s.equals(caller.getSymbol()))
			return caller.getName();
			
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
					VCFHeaderLine vcfFormatHeaderLine=new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String,"Genotype" );
					mInfoMetaData.add(vcfFormatHeaderLine);
					vcfFormatHeaderLine=new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer,"Total Depth" );
					mInfoMetaData.add(vcfFormatHeaderLine);
					VCFFilterHeaderLine vcfFilterHeaderLine=new VCFFilterHeaderLine("LowQual","Low confidence call");
					mInfoMetaData.add(vcfFilterHeaderLine);
					vcfFilterHeaderLine=new VCFFilterHeaderLine("PASS","high confidence call");
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
