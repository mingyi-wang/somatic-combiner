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
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.*;
public class SomaticCombiner {

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
		
		CommandLineParser parser=new DefaultParser();
		Options options=new Options();
		options.addOption("l", "lofreq-snv", true, "Lofreq SNV VCF file");
		options.addOption("o", "output", true, "Output VCF file");
		options.addOption("s","strelka-snv",true,"Strelka SNV VCF file");
		options.addOption("S","strelka-indel",true,"Strelka INDEL VCF file");
		options.addOption("u","muse",true,"Muse VCF file");
		options.addOption("m","mutect",true,"Mutect VCF file");
		options.addOption("M","mutect2",true,"Mutect2 VCF file");
		options.addOption("L","Lofreq-indel",true,"Lofreq INDEL VCF file");
		options.addOption("v","vardict",true,"Vardict VCF file");
		String vardictFilePath=null;
		String outputFilePath=null;
		
		String museFilePath=null;
		CommandLine line = null;
		try{
			line=parser.parse(options, args);
		}
		catch (ParseException e) {
            System.out.println(e.getMessage());
            //formatter.printHelp("utility-name", options);

            System.exit(1);
        }
		if (line.hasOption("v")) {
			vardictFilePath=line.getOptionValue("v");
		}
		if (line.hasOption("u")) {
			museFilePath=line.getOptionValue("u");
		}
		List<Variant> list =	new ArrayList<Variant>();
		String lofreqSNVFilePath=line.getOptionValue("l");
		String lofreqIndelFilePath=line.getOptionValue("L");
		String strelkaSNVFilePath=line.getOptionValue("s");
		String strelkaIndelFilePath=line.getOptionValue("S");
		
		museFilePath=line.getOptionValue("u");
		String mutectFilePath=line.getOptionValue("m");
		String mutect2FilePath=line.getOptionValue("M");
		vardictFilePath=line.getOptionValue("v");
		outputFilePath=line.getOptionValue("o");
		VCFFile vardictVCF=new VCFFile(vardictFilePath,"vardict");
		
		VCFFile museVCF=new VCFFile(museFilePath,"muse");
		VCFFile mutectVCF=new VCFFile(mutectFilePath,"mutect");
		VCFFile mutect2VCF=new VCFFile(mutect2FilePath,"mutect2");
		VCFFile strelkaSNV=new VCFFile(strelkaSNVFilePath,"strelka");
		VCFFile strelkaINDEL=new VCFFile(strelkaIndelFilePath,"strelka");
		VCFFile lofreqSNV=new VCFFile(lofreqSNVFilePath,"lofreq");
		VCFFile lofreqINDEL=new VCFFile(lofreqIndelFilePath,"lofreq");
		museVCF.importVariants(list);
		vardictVCF.importVariants(list);
		mutectVCF.importVariants(list);
		mutect2VCF.importVariants(list);
		strelkaSNV.importVariants(list);
		strelkaINDEL.importVariants(list);
		lofreqSNV.importVariants(list);
		lofreqINDEL.importVariants(list);
		System.out.println(vardictFilePath);
		List<Variant> mergedList=new ArrayList<Variant>();
		System.out.println("Starting merging...");
		long i=0;
		for (Variant p : list) {
			i++;
			System.out.println(i);
			if (i==10475)
				System.out.println("found!");
			int index=mergedList.indexOf(p);
			if (index!=-1) {
				Variant mergedVariant=mergedList.get(index).merge(p);
				if (mergedVariant==null)
					mergedList.add(p);
				else
				    mergedList.set(index, mergedVariant);
			}
			else {
				mergedList.add(p);
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
		museVCF.mergeHeader(vardictVCF.getmHeader());
		museVCF.mergeHeader(mutectVCF.getmHeader());
		museVCF.mergeHeader(mutect2VCF.getmHeader());
		museVCF.mergeHeader(strelkaSNV.getmHeader());
		museVCF.mergeHeader(strelkaINDEL.getmHeader());
		museVCF.mergeHeader(lofreqSNV.getmHeader());
		museVCF.mergeHeader(lofreqINDEL.getmHeader());
		VCFFile.writeHeader(museVCF.getmHeader(), bw);
		i=0;
		for (Variant p : mergedList) {
			i++;
			System.out.println("Writing:"+i);
			if (i==45)
				System.out.println("found!");
			museVCF.writeVariants(bw,p);
		}
		bw.close();
		
	}

}
