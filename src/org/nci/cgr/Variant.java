package org.nci.cgr;


import java.util.ArrayList;
import java.util.HashMap;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.Comparator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class Variant {
        private VariantContext variantContext;
	    // private String variantKey;
	    private String caller;
		private Byte set;
		private int priority;
	    public Variant( VariantContext variant, String callerID) throws ClassNotFoundException{
	    	variantContext=variant;    	
	    	caller=callerID;
	    	switch(caller) {
	    	  case "vardict":
	    	    set=0b000001;
	    	    priority=3;
	    	    break;
	    	  case "muse":
	    		set=0b000010;
	    		priority=4;
	    		break;
	    	  case "strelka":
	    		set=0b000100;
	    		priority=2;
	    		break;
	    	  case "lofreq":
	    		set=0b001000;
	    		priority=1;
	    		break;
	    	  case "mutect2":
	    		set=0b010000;
	    		priority=6;
	    		break;
	    	  case "mutect":
	    		set=0b100000;
	    		priority=5;
	    		break;
	    	}
	    	
	    	// variantKey=variantID;
	    	
//	    	String myDriver = "org.gjt.mm.mysql.Driver";
//			String myUrl = "jdbc:mysql://localhost/annotation_db";
//			Class.forName(myDriver);
//			connection = DriverManager.getConnection(myUrl, "root", "zh990122");
			
	    }
	    public Variant( VariantContext variant, String callerID,int p,Byte s) throws ClassNotFoundException{
	    	variantContext=variant;    	
	    	caller=callerID;
	    	priority=p;
	    	set=s;
	    }
	    	
	    public static Comparator<Variant> VariantComparator = new Comparator<Variant>() {
	    	public int compare(Variant v1, Variant v2 ) {
	    	   String contig1 = v1.getVariantContext().getContig();
	    	   String contig2 = v2.getVariantContext().getContig();
	    
	    	   if (contig1.startsWith("chr"))
	    		   contig1=contig1.substring(3);
	    	   if (contig2.startsWith("chr"))
	    		   contig2=contig2.substring(3);
	    	   
	    	   //ascending order
	    	   if (contig1.matches("\\d+") && contig2.matches("\\d+")) {
	    		   int contigNum1=Integer.parseInt(contig1);
		    	   int contigNum2=Integer.parseInt(contig2);
	    		   if (contigNum1==contigNum2)
	    			   return Integer.compare(v1.getVariantContext().getStart(),v2.getVariantContext().getStart());
	    		   else
	    			   return Integer.compare(contigNum1,contigNum2);
	    	   }
	    	   else {
	    		   if(contig1.equals(contig2))
	    		      return Integer.compare(v1.getVariantContext().getStart(),v2.getVariantContext().getStart());
	    		   else
	    			  return contig1.compareTo(contig2);
	    			   
	    	   }
	    	   
	        }};
	    
	    public VariantContext getVariantContext() {
			return variantContext;
		}

		public void setVariantContext(VariantContext variantContext) {
			this.variantContext = variantContext;
		}

		
		public String getCaller() {
			return caller;
		}
		
		

		public Byte getSet() {
			return set;
		}
		
		public void setSet(Byte set) {
			this.set = set;
		}
		public void setCaller(String caller) {
			this.caller = caller;
		}

		@Override
	    public int hashCode() {
			return Objects.hash(variantContext.getContig(), variantContext.getStart());
	    		
	    }
	    
		@Override
		public boolean equals(Object o) {
	        if (this == o) return true;
	        if (o == null || getClass() != o.getClass()) return false;
	        Variant that = (Variant) o;
	        return this.variantContext.getContig().equals(that.variantContext.getContig()) &&
	                this.variantContext.getStart()==that.variantContext.getStart();
	    }

	    public Variant merge(Variant other) throws ClassNotFoundException {
	    	assert(this.equals(other));
	    	VariantContextBuilder build=new VariantContextBuilder();
	    	List<Allele> alleles = this.variantContext.getAlleles();
	    	Boolean multipleCallerCalled=false;
	    	if (this.getCaller().contains(","))
	    		multipleCallerCalled=true;
			if (alleles.equals(other.variantContext.getAlleles())) {
				Byte s=(byte) (this.set|other.set);
				Map<String,Object> thisInfo=new HashMap<String,Object>();
				Map<String,Object> otherInfo=new HashMap<String,Object>();
				Map<String,Object> mergedInfo=new HashMap<String,Object>();
				
				thisInfo=this.variantContext.getAttributes();

				for (String key:thisInfo.keySet()){
					if (multipleCallerCalled)
						mergedInfo.put(key, thisInfo.get(key));
					else
					    mergedInfo.put(this.caller+"_"+key, thisInfo.get(key));
				}
				otherInfo=other.variantContext.getAttributes();
				for (String key:otherInfo.keySet())
					mergedInfo.put(other.caller+"_"+key, otherInfo.get(key));
				
								
				GenotypesContext thisGenotypes=null;
				GenotypesContext otherGenotypes=null;
				List<Genotype> mergedGenotypes=new ArrayList<Genotype>();
				if (this.variantContext.hasGenotypes() && other.variantContext.hasGenotypes()) {
				   thisGenotypes=this.variantContext.getGenotypes();
				   otherGenotypes=other.variantContext.getGenotypes();
				   if (thisGenotypes.size()==otherGenotypes.size()) {
						Set<String> thisSampleNames = thisGenotypes.getSampleNames();
						Set<String> otherSampleNames = otherGenotypes.getSampleNames();
						if (thisSampleNames.size() == 2) {
							String[] thisSamples = new String[2];
							boolean thisFoundTumor = false;
							for (String sampleName : thisSampleNames)
								if (sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG1)||sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG2)) {
									thisSamples[0] = sampleName;
									thisFoundTumor = true;
								} else
									thisSamples[1] = sampleName;
							String[] otherSamples = new String[2];
							boolean otherFoundTumor = false;
							for (String sampleName : otherSampleNames)
								if (sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG1)||sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG2)) {
									otherSamples[0] = sampleName;
									otherFoundTumor = true;
								} else
									otherSamples[1] = sampleName;
	                        if (thisFoundTumor && otherFoundTumor) {
								for (int i = 0; i < thisGenotypes.size(); i++) {
									String sampleName=VCFFile.TUMOR_TAG1;
									if (i==1) sampleName=VCFFile.NORMAL_TAG;
									GenotypeBuilder mergedGenotypeBuilder = new GenotypeBuilder();
									
									Map<String, Object> mergedExtendedGT = new HashMap<String, Object>();
									Genotype tmpGenotype = null;
									Map<String, Object> thisExtendedGT = thisGenotypes.get(thisSamples[i]).getExtendedAttributes();
									Map<String, Object> otherExtendedGT = otherGenotypes.get(otherSamples[i]).getExtendedAttributes();
									for (String key : thisExtendedGT.keySet())
										if (multipleCallerCalled)
											mergedExtendedGT.put(key, thisExtendedGT.get(key));
										else
											mergedExtendedGT.put(this.caller + "_" + key, thisExtendedGT.get(key));
									for (String key : otherExtendedGT.keySet())
										mergedExtendedGT.put(other.caller + "_" + key, otherExtendedGT.get(key));
									mergedGenotypeBuilder.attributes(mergedExtendedGT);
		
									if (this.priority > other.priority)
										tmpGenotype = thisGenotypes.get(thisSamples[i]);
									else
										tmpGenotype = otherGenotypes.get(otherSamples[i]);
									mergedGenotypeBuilder.alleles(tmpGenotype.getAlleles());
									if (thisGenotypes.get(thisSamples[i]).hasAD() && otherGenotypes.get(otherSamples[i]).hasAD())
										mergedGenotypeBuilder.AD(tmpGenotype.getAD());
									else {
										if (thisGenotypes.get(thisSamples[i]).hasAD())
											mergedGenotypeBuilder.AD(thisGenotypes.get(thisSamples[i]).getAD());
										if (otherGenotypes.get(otherSamples[i]).hasAD())
											mergedGenotypeBuilder.AD(otherGenotypes.get(otherSamples[i]).getAD());
									}
		
									if (thisGenotypes.get(thisSamples[i]).hasDP() && otherGenotypes.get(otherSamples[i]).hasDP())
										mergedGenotypeBuilder.DP(tmpGenotype.getDP());
		
									else {
										if (thisGenotypes.get(thisSamples[i]).hasDP())
											mergedGenotypeBuilder.DP(thisGenotypes.get(thisSamples[i]).getDP());
										if (otherGenotypes.get(i).hasDP())
											mergedGenotypeBuilder.DP(otherGenotypes.get(otherSamples[i]).getDP());
									}
		                            mergedGenotypeBuilder.name(sampleName);
									mergedGenotypes.add(mergedGenotypeBuilder.make());
								}
						    }
						}
				   }
				   else {
					   List<Genotype> tmpGenotypes=null;
						if (this.variantContext.getGenotypes().size()>other.variantContext.getGenotypes().size()) 
							tmpGenotypes=this.variantContext.getGenotypes();
						else
							tmpGenotypes=other.variantContext.getGenotypes();
						if (tmpGenotypes!=null) {
							for (Genotype tmpGenotype: tmpGenotypes)
						       mergedGenotypes.add(tmpGenotype);
						}
				   }
				   
				}else {
					List<Genotype> tmpGenotypes=null;
					if (this.variantContext.hasGenotypes()) 
						tmpGenotypes=this.variantContext.getGenotypes();
					
					if (other.variantContext.hasGenotypes()) 
						tmpGenotypes=other.variantContext.getGenotypes();
					if (tmpGenotypes!=null) {
						for (Genotype tmpGenotype: tmpGenotypes)
					       mergedGenotypes.add(tmpGenotype);
					}
				}
				
				
				build.alleles(this.variantContext.getAlleles());
				build.chr(this.variantContext.getContig());
				build.start(this.variantContext.getStart());
				if (this.variantContext.getID().equals(other.variantContext.getID()))
				   build.id(this.variantContext.getID());
				else {
				   if (this.variantContext.getID().equals("."))
					   build.id(other.variantContext.getID());
				   if (other.variantContext.getID().equals("."))
					   build.id(this.variantContext.getID());
				}   
				build.stop(this.variantContext.getEnd());
				build.genotypes(mergedGenotypes);
				build.attributes(mergedInfo);
				VariantContext mergedVc=build.make();
				int mergedPriority=this.priority;
				if (this.priority<other.priority)
					mergedPriority=other.priority;
				return new Variant(mergedVc,this.caller+","+other.caller,mergedPriority,s);
			  
			}
			return null;
	    	
	    }
	 }
