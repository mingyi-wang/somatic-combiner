package org.nci.cgr;

import java.io.BufferedWriter;
import java.io.File;

import java.io.IOException;
import java.io.Writer;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;



import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.*;

public class VCFFile {
	private String filePath;
	private String caller;
	private File vcfFile;
	private static VCFFileReader vcfFileReader;
	private VCFHeader mHeader = null;
	public static String TUMOR_TAG1="TUMOR";
	public static String TUMOR_TAG2="TUMOUR";
	public static final String NORMAL_TAG="NORMAL";
	private static final String VERSION_LINE =
            VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF4_2.getFormatString() + "=" + VCFHeaderVersion.VCF4_2.getVersionString();
	public VCFFile(String filePath, String caller) {
		super();
		this.filePath = filePath;
		this.caller = caller;
		
	}
	
	
	public VCFHeader getmHeader() {
		return mHeader;
	}


	public void setmHeader(VCFHeader mHeader) {
		this.mHeader = mHeader;
	}


	public void writeVariants(BufferedWriter bw,Variant vv, int snvCallerNum,int indelCallerNum) throws IOException {
//		if (vv.getVariantContext().getStart()==7882066)
//			System.out.println("found!");
		bw.append(vv.getVariantContext().getContig()+"\t"+vv.getVariantContext().getStart()+"\t"+vv.getVariantContext().getID()+"\t"+
	              vv.getVariantContext().getAlleles().get(0).getBaseString()+"\t");
		
		if (vv.getVariantContext().getAlleles().size()>1)
			bw.append(vv.getVariantContext().getAlleles().get(1).getBaseString()+"\t");
		else
			bw.append(vv.getVariantContext().getAlleles().get(0).getBaseString()+"\t");
		
		Map<String,Object> info=vv.getVariantContext().getAttributes();
		String infoContent="";
		String qualContent=".";
		for (String key:info.keySet()) {
			if (key=="qual")
				qualContent=info.get(key).toString();
			else
			    infoContent=infoContent+key+"="+info.get(key).toString()+";";
		}
		bw.append(qualContent+"\t");
		
		if (vv.getVariantContext().isSNP()) {
			if ((float)Integer.bitCount(vv.getSet()) / snvCallerNum >= 0.5)
				bw.append("PASS" + "\t");
			else
				bw.append("LowQual" + "\t");
		} else {
			if ((float)Integer.bitCount(vv.getSet()) / indelCallerNum >= 0.5)
				bw.append("PASS" + "\t");
			else
				bw.append("LowQual" + "\t");
		}
						
		bw.append(infoContent+"NumTools="+Integer.bitCount(vv.getSet())+";"+String.format("%6s", Integer.toBinaryString(vv.getSet() & 0xFF)).replace(' ', '0'));
		String format="";
		String gtString="";
		//vcfEncoder.write(bw, vv.getVariantContext());
		if (vv.getVariantContext().hasGenotypes()) {
		   Set<String> sampleNames=vv.getVariantContext().getSampleNames();
		   if ((sampleNames.size()!=2) && (sampleNames.size()!=0)) {
			   System.err.println("Error: the samples in "+filePath+" not equals to 2 or 0! The variants in this file will be skipped!");
               return;
		   }
		   if (sampleNames.size()==2) {
			   String[] samples=new String[2];
			   boolean foundTumor=false;
			   for(String sampleName:sampleNames) 
				   if (sampleName.toUpperCase().contains("TUMOR")||sampleName.toUpperCase().contains("TUMOUR")) {
					  samples[0]=sampleName;
					  foundTumor=true;  
				   }
				   else
					  samples[1]=sampleName;
			   if (!foundTumor) {
				   System.err.println("Error: the samples in "+filePath+" does not contain TUMOR in the genotype columns! The variants in this file will be skipped!");
	               return; 
			   }
			   for(int i=0;i<samples.length;i++) {
				   Genotype gt=vv.getVariantContext().getGenotype(samples[i]);
				   
				   if (i==0) 
					   if (gt.getAlleles().size()>0) format="GT";
				   if (gt.getAlleles().size()==1)
					   gtString=gtString+gt.getAlleles().get(0).getBaseString()+"/"+gt.getAlleles().get(0).getBaseString();
				   else
					   if (gt.getAlleles().size()==2)
				         gtString=gtString+gt.getAlleles().get(0).getBaseString()+"/"+gt.getAlleles().get(1).getBaseString();
				   if (gt.hasDP()) {
					   if (i==0) {
						   if (format.length()==0)
							   format="DP";
						   else
							   format=format+":DP";
					   }
					   if (gt.getAlleles().size()==0)
						   gtString=gtString+gt.getDP();
					   else
				           gtString=gtString+":"+gt.getDP();
				   }
			
				   if (gt.hasAD()) {  
					   gtString=gtString+":"+Arrays.toString(gt.getAD()).replace("[", "").replace("]", "").replace(" ","").trim();
					   if (i==0) format=format+":AD";   
				   }
				   Map<String,Object> extendedAttributes=gt.getExtendedAttributes() ;
				   for (String key:extendedAttributes.keySet()) {
					   if (i==0) format=format+":"+key;
					   gtString=gtString+":"+extendedAttributes.get(key).toString();
				   }
				   
				   gtString=gtString+"\t";
			   }   
			}
		}
		bw.append("\t"+format);
		while (gtString.endsWith("\t")) {
			gtString = gtString.substring(0, gtString.length()-1);
			bw.append("\t"+gtString);
		}
		
		bw.newLine();
		bw.flush();
	}
	
	private static void rejectVCFV43Headers(final VCFHeader targetHeader) {
        if (targetHeader.getVCFHeaderVersion() != null && targetHeader.getVCFHeaderVersion().isAtLeastAsRecentAs(VCFHeaderVersion.VCF4_3)) {
            throw new IllegalArgumentException(String.format("Writing VCF version %s is not implemented", targetHeader.getVCFHeaderVersion()));
        }

    }
	
	public void mergeHeader(VCFHeader header) {
		for (final VCFHeaderLine line : header.getMetaDataInSortedOrder()) {
			if (VCFHeaderVersion.isFormatString(line.getKey()))
				continue;

			mHeader.addMetaDataLine(line);
		}
		
	}
	
	public static VCFHeader writeHeader(VCFHeader header, final Writer writer) {

		try {
			rejectVCFV43Headers(header);
			String versionLine=VERSION_LINE;

			// the file format field needs to be written first
			writer.write(versionLine + "\n");

			for (final VCFHeaderLine line : header.getMetaDataInSortedOrder()) {
				if (VCFHeaderVersion.isFormatString(line.getKey()))
					continue;

				writer.write(VCFHeader.METADATA_INDICATOR);
				writer.write(line.toString());
				writer.write("\n");
			}

			// write out the column line
			writer.write(VCFHeader.HEADER_INDICATOR);
			boolean isFirst = true;
			for (final VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) {
				if (isFirst)
					isFirst = false; // don't write out a field separator
				else
					writer.write(VCFConstants.FIELD_SEPARATOR);
				writer.write(field.toString());
			}

			if (header.hasGenotypingData()) {
				writer.write(VCFConstants.FIELD_SEPARATOR);
				writer.write("FORMAT");
				for (final String sample : header.getGenotypeSamples()) {
					writer.write(VCFConstants.FIELD_SEPARATOR);
					writer.write(sample);
				}
			}

			writer.write("\n");
			writer.flush(); // necessary so that writing to an output stream will work
		} catch (IOException e) {
			throw new RuntimeIOException("IOException writing the VCF header to " , e);
		}

		return header;
	}
	
	public int importVariants(List<Variant> list) throws ClassNotFoundException, SQLException {
		vcfFile = new File(filePath);

		vcfFileReader = new VCFFileReader(vcfFile, false);
		
		mHeader = vcfFileReader.getFileHeader();
		ArrayList<String> sampleNames=mHeader.getSampleNamesInOrder();
		boolean foundTumor=false;
		if (sampleNames.size()==2) {
			for(String s:sampleNames) {
				if (s.toUpperCase().contains(TUMOR_TAG1) || s.toUpperCase().contains(TUMOR_TAG2))
				  foundTumor=true;
			}
			if (!foundTumor) {
				System.out.println("Warning: "+filePath+" is skipped due to no TUMOR or TUMOUR found in the samplename line!");
				return 0;
			}
		}
		else
			if (sampleNames.size()>2) {
				System.out.println("Warning: "+filePath+" is skipped due to more than two samples in the VCF!");
				return 0;
			}
				
		// System.out.println(ss);
		final CloseableIterator<VariantContext> variantIterator = vcfFileReader.iterator();
		while (variantIterator.hasNext()) {
			final VariantContext vc = variantIterator.next();
//			if (vc.getStart()==61851)
//				System.out.println("found!");
			System.out.println(caller+" "+vc.getContig()+":"+vc.getStart());
			// Set<String> filters=vc.getFilters();
			Boolean ff = vc.isFiltered();
			Set<String> filters = vc.getFilters();
			String filterContent = "";

			if (!filters.isEmpty()) {
				Iterator<String> iter = filters.iterator();
				while (iter.hasNext()) {
					filterContent = filterContent + caller + "_" + iter.next() + ",";
				}
				filterContent = filterContent.substring(0, filterContent.length() - 1);
			}

			// Set<String> f2=kv.getFiltersMaybeNull();
			//
			// System.out.println(ff);
			// Boolean passed=false;
			// for (String f:filters) {
			// System.out.println(f);
			// if (f.equals("PASS")) {
			// passed=true;
			// break;
			// }
			// }

			if (ff == false) {
				if (vc.isBiallelic()) {
					Variant variant=new Variant(vc,caller);
					list.add(variant);
				}
				else {
					//Muse split allele and genotypes
					List<Allele> alleles = vc.getAlleles();
					
					int alleleCount=alleles.size();
					Variant variant=new Variant(vc,caller);
					for (int i=1;i<alleleCount;i++) {
					  VariantContextBuilder build=new VariantContextBuilder();
					  VariantContext splitVC=null;
					  VariantContext tmpVC=variant.getVariantContext();
					  List<Allele> tmpAlleles=tmpVC.getAlleles();
					  List<Genotype> tmpGenotypes=tmpVC.getGenotypes();
					 
					  for (int j=1;j<alleleCount;j++) 
						 if ((j!=i) && (j!=0)) 
					       tmpAlleles.remove(j);
					  for (int j=0;j<tmpAlleles.size();j++) 
						  System.out.println(tmpAlleles.get(j).getBaseString());
						  				  
					  List<Genotype> splitGenotypes=new ArrayList<Genotype>();
					  for (Genotype gt:tmpGenotypes) {
						  if (gt.hasAD()) {
							  GenotypeBuilder splitGenotypeBuilder=new GenotypeBuilder();
							  int[] ADs=gt.getAD();
							  int[] newADs=new int[2];
							  for (int j=0,k=0;j<alleleCount;j++) 
									 if ((j==i) || (j==0)) 
								       newADs[k++]=ADs[j];									   
							  
							  List<Allele> tmpGtAlleles=gt.getAlleles();
							  
							  for (int j=0;j<tmpGtAlleles.size();j++) {
								  boolean found=false;
								  for (int k=0;k<tmpAlleles.size();k++) {
									  if (tmpGtAlleles.get(j).equals(tmpAlleles.get(k))) {
										  found=true;
										  break;
									  }
								  }
								  if (!found) {
									  tmpGtAlleles.remove(j); 
									  break;
								  }  
							  }
							  for (int j=0;j<tmpGtAlleles.size();j++) 
								  System.out.println("genotype:"+tmpGtAlleles.get(j).getBaseString());
								  							  
							  splitGenotypeBuilder.alleles(tmpGtAlleles);
							  splitGenotypeBuilder.AD(newADs);
							  if (gt.hasDP()) splitGenotypeBuilder.DP(gt.getDP());							 
							  Map<String,Object> gtExtendedAttributes=gt.getExtendedAttributes();
						      for(String key:gtExtendedAttributes.keySet()) {
						    	  String cName=gtExtendedAttributes.get(key).getClass().getName();
//								  System.out.print(cName);
								  if (cName.contains("ArrayList")){
										ArrayList fieldsForKey=(ArrayList) gtExtendedAttributes.get(key);
										if (fieldsForKey.size()==alleleCount) 
											for (int j=1;j<alleleCount;j++) 
												 if ((j!=i) && (j!=0)) 
													 fieldsForKey.remove(j);																		
								  }
						      }
						      splitGenotypeBuilder.name(gt.getSampleName());
						      splitGenotypeBuilder.attributes(gtExtendedAttributes);
						      splitGenotypes.add(splitGenotypeBuilder.make());
						  }
						  else 
							  splitGenotypes.add(gt);
						  						 						 
					  }
					  build.alleles(tmpAlleles);					  
					  build.chr(vc.getContig());
					  build.start(vc.getStart());
					  build.id(vc.getID());
					  if (vc.hasLog10PError()) 
						  build.log10PError(vc.getLog10PError());  
					  build.stop(vc.getEnd());
					  build.genotypes(splitGenotypes);
					  build.attributes(vc.getAttributes());
					  splitVC=build.make();
					  Variant splitVariant=new Variant(splitVC,caller);
					  list.add(splitVariant);
					}
				}
			    	
			}
    	}
		return 1;
	}
	
	
}
