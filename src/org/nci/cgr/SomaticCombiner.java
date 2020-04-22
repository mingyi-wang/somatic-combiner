package org.nci.cgr;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.logging.StreamHandler;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
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
	public static Logger logger = Logger.getLogger(SomaticCombiner.class.getName());
	public static String version="V1.02";
	public static ArrayList<String> medianAFs=new ArrayList<String>();
	public static void main(String[] args) throws ClassNotFoundException, SQLException, IOException {
		InputStream is=SomaticCombiner.class.getClassLoader().getResourceAsStream("mylogging.properties");
		LogManager.getLogManager().readConfiguration(is);
		logger.setLevel(Level.FINE);
        // logger.addHandler(new ConsoleHandler());
        // logger.addHandler(new StreamHandler(System.out,new MyFormatter()));
		initCallers();
		CommandLineParser parser=new DefaultParser();
		Options options=new Options();
		for (Caller caller :callerList)
			options.addOption(caller.getSymbol(),caller.getName(),true,caller.getDescription());
		
		options.addOption("o","output",true,"Output VCF file");
		options.addOption("h","help",false,"print this message");
		String outputFilePath=null;
		

		CommandLine line = null;
		try{
			line=parser.parse(options, args);
		}
		catch (ParseException e) {
			logger.log(Level.SEVERE, e.getMessage());
            //formatter.printHelp("utility-name", options);

            System.exit(1);
        }
		int count=0;
		if (line.hasOption("h")) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "somaticCombiner "+version, options );
			System.exit(1);	
		}
		for (Caller caller:callerList) {
			if (line.hasOption(caller.getSymbol())) {
				caller.setFilePath(line.getOptionValue(caller.getSymbol()));
				File tmpFile=new File(caller.getFilePath());
				if (tmpFile.exists()) 
					count++;
				else {
					// System.out.println("Error:"+caller.getFilePath()+" does not existed! Please check for that! The program will skip this VCF file!");
					logger.log(Level.SEVERE, "Error:"+caller.getFilePath()+" does not existed! Please check for that! The program will skip this VCF file!");
					System.exit(1);
				}				
			}
		}
		if (count<2) {
			if (count==1)
			    logger.log(Level.SEVERE,"Error: Only one VCF is available. No merge needed!");
			else {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "somaticCombiner "+version, options );
			}
			System.exit(1);
		}
		outputFilePath=line.getOptionValue("o");
		
		int snvCallerNum=0;
		int indelCallerNum=0;
		List<Variant> list =new ArrayList<Variant>();
		logger.log(Level.INFO,"somaticCombiner "+version);
		System.out.println("somaticCombiner "+version);
		Option[] allOptions=line.getOptions();
		String parameters="";
		for (Option o:allOptions) {
			System.out.print("-"+o.getOpt()+" "+o.getValue()+" ");
			parameters+="-"+o.getOpt()+" "+o.getValue()+" ";
		}
		System.out.println();
		logger.log(Level.INFO,"Parameters: "+parameters);
		logger.log(Level.INFO,"Loading VCFs...");
		for (Caller caller:callerList) {
			if (caller.getFilePath()!=null) {
				VCFFile vcfFile=new VCFFile(caller.getFilePath(),caller.getName(),caller.getPriority(),caller.getSet());
				logger.log(Level.INFO,"Loading "+caller.getName()+" "+caller.getFilePath()+" ...");
//				if (caller.getName().contains("strelka-indel")) {
//					System.out.println("Found!");
//				}
				int cnt=vcfFile.importVariants(list);
			    if(cnt>=0) {
			    	logger.log(Level.INFO,"Finished loading "+caller.getName()+" from "+caller.getFilePath()+". Total count:"+cnt+".");
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
			    else 
			    	logger.log(Level.SEVERE,"Warning: Loading "+caller.getName()+" "+caller.getFilePath()+" is failed!");
			    	
			    
			    
			}
		}
		if (callerSymbols.length()<=1) {
			logger.log(Level.SEVERE,"Caller number is less than 2. Exit!");
            System.exit(1);
		}
		callerNames=callerNames.substring(0, callerNames.length()-1);
		List<MergedVariant> mergedList=new ArrayList<MergedVariant>();
		logger.log(Level.INFO,"Starting merging...");
		long i=0;
		for (Variant p : list) {
			i++;
			if (i%10000==0)
			  logger.log(Level.INFO,"Processed: "+i);
//			if (p.getVariantContext().getStart()==30880432)
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
		logger.log(Level.INFO,"Finished merging VCFs. Processed: "+i+" in total.");
		logger.log(Level.INFO,"Sorting merged variants ...");
		Collections.sort(mergedList,Variant.VariantComparator);
		logger.log(Level.INFO,"Writing VCF ...");

		
		BufferedWriter bw=new BufferedWriter(new FileWriter(new File(outputFilePath)));
		VCFFile firstVCFFile=prepareHeader();	
		VCFFile.writeHeader(firstVCFFile.getmHeader(), bw);
				
		i=0;
		for (Variant p : mergedList) {
			i++;
//			if (p.getVariantContext().getStart()==1717242)
//				System.out.println("Output");
			if (i%10000==0)
		    	logger.log(Level.INFO,"Finished writing "+i);
			firstVCFFile.writeVariants(bw,p,snvCallerNum,indelCallerNum);
		}
		logger.log(Level.INFO,"Finished writing the merged VCF and processed "+i+" in total.");
		bw.close();
		logger.log(Level.INFO,"Writing VCF is done!");
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
