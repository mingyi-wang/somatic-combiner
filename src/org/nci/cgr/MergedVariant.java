package org.nci.cgr;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;


import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class MergedVariant extends Variant   
{
	 public MergedVariant(VariantContext variant, String callerID,int p,Byte s) throws ClassNotFoundException{

	    	super(variant,callerID,p,s);
	    }
	@Override
    public int hashCode() {
		String alt="";
		if (variantContext.getAlleles().size()==2)
			alt=variantContext.getAlleles().get(1).getBaseString();
		else
			alt=variantContext.getAlleles().get(0).getBaseString();
		return Objects.hash(variantContext.getContig(), variantContext.getStart(),variantContext.getReference().getBaseString(),alt);
    		
    }
    
	@Override
	public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Variant that = (Variant) o;
        return this.variantContext.getContig().equals(that.variantContext.getContig()) &&
                this.variantContext.getStart()==that.variantContext.getStart() && this.variantContext.getAlleles().equals(that.variantContext.getAlleles());
    }

    public MergedVariant merge(MergedVariant other) throws ClassNotFoundException {
    	assert(this.equals(other));
    	VariantContextBuilder build=new VariantContextBuilder();
    	// List<Allele> alleles = this.variantContext.getAlleles();
    	Boolean multipleCallerCalled=false;
    	if (this.getCaller().contains(",")) multipleCallerCalled=true;
		//if (alleles.equals(other.variantContext.getAlleles())) {
		Byte s=(byte) (this.set|other.set);
		if (s==this.set) return(null);  //duplicate variant from a same caller
		Map<String,Object> thisInfo=new HashMap<String,Object>();
		Map<String,Object> otherInfo=new HashMap<String,Object>();
		Map<String,Object> mergedInfo=new HashMap<String,Object>();
		
		thisInfo=this.variantContext.getAttributes();

		for (String key:thisInfo.keySet()){
			if (multipleCallerCalled)
				mergedInfo.put(key, thisInfo.get(key));
			else
			    mergedInfo.put(SomaticCombiner.callerName(this.caller)+"_"+key, thisInfo.get(key));
		}
		otherInfo=other.variantContext.getAttributes();
		for (String key:otherInfo.keySet())
			mergedInfo.put(SomaticCombiner.callerName(other.caller)+"_"+key, otherInfo.get(key));
		
		String qual="";
		if (this.getVariantContext().hasLog10PError()) {
			if (!multipleCallerCalled) {
			qual=String.valueOf(this.getVariantContext().getPhredScaledQual());
			mergedInfo.put(SomaticCombiner.callerName(this.caller)+"_"+VCFFile.QUAL_TAG, qual);
			}
		}
		if (other.getVariantContext().hasLog10PError()) {
			qual=String.valueOf(other.getVariantContext().getPhredScaledQual());
			mergedInfo.put(SomaticCombiner.callerName(other.caller)+"_"+VCFFile.QUAL_TAG, qual);
		}
		
						
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
							for (String key : thisExtendedGT.keySet()) {
//								if (key.contains("FDP")) 
//									System.out.println("Found!");
								
								if (multipleCallerCalled)
									mergedExtendedGT.put(key, thisExtendedGT.get(key));
								else
									mergedExtendedGT.put(SomaticCombiner.callerName(this.caller) + "_" + key, thisExtendedGT.get(key));
							}
							for (String key : otherExtendedGT.keySet())
								mergedExtendedGT.put(SomaticCombiner.callerName(other.caller) + "_" + key, otherExtendedGT.get(key));
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
					for (Genotype tmpGenotype:tmpGenotypes)
					
				       mergedGenotypes.add(tmpGenotype);
				}
		   }
		   
		}else {
			GenotypesContext tmpGenotypes=null;
			Set<String> tmpSampleNames = null;
			String tmpCaller="";
			if (this.variantContext.hasGenotypes()) { 
				tmpGenotypes=this.variantContext.getGenotypes();
			    tmpSampleNames=this.variantContext.getSampleNames();
			    tmpCaller=this.getCaller();
			}
			if (other.variantContext.hasGenotypes()) {
				tmpGenotypes=other.variantContext.getGenotypes();
			    tmpSampleNames=other.variantContext.getSampleNames();
			    tmpCaller=other.getCaller();
			}
			
			if (tmpGenotypes!=null) {
				Map<String, Object> mergedExtendedGT = new HashMap<String, Object>();
				if (tmpSampleNames.size() == 2) {
					String[] tmpSamples = new String[2];
					boolean tmpFoundTumor = false;
					for (String sampleName : tmpSampleNames)
						if (sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG1)||sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG2)) {
							tmpSamples[0] = sampleName;
							tmpFoundTumor = true;
						} else
							tmpSamples[1] = sampleName;
					if (tmpFoundTumor) {
						
					    for (int i = 0; i < tmpGenotypes.size(); i++) {
					    	GenotypeBuilder mergedGenotypeBuilder = new GenotypeBuilder();
					    	String sampleName=VCFFile.TUMOR_TAG1;
							if (i==1) sampleName=VCFFile.NORMAL_TAG;
				    	    Map<String, Object> thisExtendedGT = tmpGenotypes.get(tmpSamples[i]).getExtendedAttributes();
					        for (String key : thisExtendedGT.keySet()) 
						    	mergedExtendedGT.put(SomaticCombiner.callerName(tmpCaller) + "_" + key, thisExtendedGT.get(key));
					    
						    mergedGenotypeBuilder.attributes(mergedExtendedGT);
	
					        mergedGenotypeBuilder.alleles(tmpGenotypes.get(tmpSamples[i]).getAlleles());
				        	if (tmpGenotypes.get(tmpSamples[i]).hasAD() )
					   	       mergedGenotypeBuilder.AD(tmpGenotypes.get(tmpSamples[i]).getAD());
				        	else {
						       if (tmpGenotypes.get(tmpSamples[i]).hasAD())
					        		mergedGenotypeBuilder.AD(tmpGenotypes.get(tmpSamples[i]).getAD());
						    }
	
				        	if (tmpGenotypes.get(tmpSamples[i]).hasDP())
						       mergedGenotypeBuilder.DP(tmpGenotypes.get(tmpSamples[i]).getDP());
			        		else {
					            if (tmpGenotypes.get(tmpSamples[i]).hasDP())
							      mergedGenotypeBuilder.DP(tmpGenotypes.get(tmpSamples[i]).getDP());
					        }
			                mergedGenotypeBuilder.name(sampleName);
					        mergedGenotypes.add(mergedGenotypeBuilder.make());
			    	
			            }
		            }
				}
			       
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
		return new MergedVariant(mergedVc,this.caller+","+other.caller,mergedPriority,s);
		  
	   //}
		//return null;
    	
    }
    
    
}
