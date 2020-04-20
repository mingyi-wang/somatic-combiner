package org.nci.cgr;


import java.util.ArrayList;
import java.util.HashMap;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.logging.Level;
import java.util.Comparator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;

public class Variant {
	    protected VariantContext variantContext;
	    // private String variantKey;
	    protected String caller;
	    protected Byte set;
	    protected int priority;
	    
	    public Variant( VariantContext variant, String callerID,int p,Byte s) throws ClassNotFoundException{
	    	variantContext=variant;    	
	    	caller=callerID;
	    	priority=p;
	    	set=s;
	    }
	    	
	    
	    public int getTumorDP() {
	    	if (variantContext.hasGenotypes()) {
	    		Set<String> sampleNames=variantContext.getSampleNames();
	    		String tumorName="";
	    		boolean foundTumor=false;
				for(String sampleName:sampleNames) 
					   if (sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG1)||sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG2)) {
						  tumorName=sampleName;
						  foundTumor=true;
						  break;
					   }
					  
				if (!foundTumor) 
		               return 0;
		        else {
		        	if (variantContext.getGenotype(tumorName).hasDP())
		        	   return variantContext.getGenotype(tumorName).getDP();
		        	else {
		        		if (variantContext.getGenotype(tumorName).hasExtendedAttribute("ALT_F1R2") && variantContext.getGenotype(tumorName).hasExtendedAttribute("ALT_F2R1")
		        				&& variantContext.getGenotype(tumorName).hasExtendedAttribute("REF_F1R2") && variantContext.getGenotype(tumorName).hasExtendedAttribute("REF_F2R1") ) {
			        		int DP=Integer.parseInt(variantContext.getGenotype(tumorName).getExtendedAttribute("ALT_F1R2").toString())+
			        		Integer.parseInt(variantContext.getGenotype(tumorName).getExtendedAttribute("ALT_F2R1").toString())+
			        		Integer.parseInt(variantContext.getGenotype(tumorName).getExtendedAttribute("REF_F1R2").toString())+
			        		Integer.parseInt(variantContext.getGenotype(tumorName).getExtendedAttribute("REF_F2R1").toString());
			        		return DP;
		        		}
		        		else {
		        			if (variantContext.getGenotype(tumorName).hasExtendedAttribute("Mutect2_ALT_F1R2") && variantContext.getGenotype(tumorName).hasExtendedAttribute("Mutect2_ALT_F2R1")
			        				&& variantContext.getGenotype(tumorName).hasExtendedAttribute("Mutect2_REF_F1R2") && variantContext.getGenotype(tumorName).hasExtendedAttribute("Mutect2_REF_F2R1") ) {
				        		int DP=Integer.parseInt(variantContext.getGenotype(tumorName).getExtendedAttribute("Mutect2_ALT_F1R2").toString())+
				        		Integer.parseInt(variantContext.getGenotype(tumorName).getExtendedAttribute("Mutect2_ALT_F2R1").toString())+
				        		Integer.parseInt(variantContext.getGenotype(tumorName).getExtendedAttribute("Mutect2_REF_F1R2").toString())+
				        		Integer.parseInt(variantContext.getGenotype(tumorName).getExtendedAttribute("Mutect2_REF_F2R1").toString());
				        		return DP;
		        			}
		        			else
		        		       return 0;
		        		
		        		}
		        		
		        	}
		        }
	    	}
	    	else {
	    		Map<String,Object> info=variantContext.getAttributes();
	    		for (String key:info.keySet()) {		
	    			if (key=="Lofreq_DP")
	    				return Integer.parseInt(info.get(key).toString());			
	    		}
	    	}
			return 0;
	    		
	    }
	    
	    public float getTumorAF() {
	    	String callerNameWithHighestPriority=SomaticCombiner.callerNameFromPriority(priority);
	    	if (callerNameWithHighestPriority.equals("Strelka") || callerNameWithHighestPriority.equals("Muse") || callerNameWithHighestPriority.equals("Varscan")){
	    		if (callerNameWithHighestPriority.equals("Muse") || callerNameWithHighestPriority.equals("Varscan") ) {
	    			if (variantContext.hasGenotype(VCFFile.TUMOR_TAG1) ) {
	    				Genotype gt=variantContext.getGenotype(VCFFile.TUMOR_TAG1);
	    				if (gt.hasAD()) {
		    			    int AD[]=gt.getAD();
		    			    if (AD.length==2)
			    			    return (float)AD[1]/ (AD[0]+AD[1]);
		    			    else {
		    			    	if (AD.length==1) {
			    			    	if (gt.hasDP()) {
			    			    		if (gt.getDP()>0)
			    			    		  return (float) AD[0]/gt.getDP();
			    			    		else {
			    			    			SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
			   			            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
			   			            		  		+ "\t"+" Varscan has DP is 0 in FORMAT!");
			   	    			    	    return -1;
			    			    		}
			    			    	}
			    			    	else {
			    			    		SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
			   			            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
			   			            		  		+ "\t"+" Varscan has no DP in FORMAT!");
			   	    			    	    return -1;
			    			    		
			    			    	}
		    			    	}
		    			    	else {
		    			    		SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
					            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
					            		  		+ "\t"+" Muse or Varscan AD lenght not equal to 2 or 1 in FORMAT!");
			    			    	return -1;
		    				    }
		    			    }
	    				}
	    				else {
	    					SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
			            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
			            		  		+ "\t"+"has no Muse AD in FORMAT!");
	    					return -1;
	    				}
	    			}
	    			else {
	    				SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
		            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
		            		  		+ "\t"+"has no Muse AD in FORMAT!");
	    				return -1;
	    			}
	    		}
	    		else {
	    			String gtRefKey=variantContext.getAlleles().get(0).getBaseString()+"U";
	    			String gtAltKey=variantContext.getAlleles().get(1).getBaseString()+"U";
	    			if (Integer.bitCount(set)>1) {
	    				gtRefKey=callerNameWithHighestPriority+"_"+gtRefKey;
	    				gtAltKey=callerNameWithHighestPriority+"_"+gtAltKey;
	    			}
	    			if (variantContext.hasGenotype(VCFFile.TUMOR_TAG1)) {
	    			   Genotype gt=variantContext.getGenotype(VCFFile.TUMOR_TAG1);
	    			   Map<String, Object> extendedGT = gt.getExtendedAttributes();
	    			   String refDP[]=null;
	    			   String altDP[]=null;
	    			   for (String key:extendedGT.keySet()) {
	            		  if (key.equals(gtRefKey)) 
                             refDP= extendedGT.get(key).toString().split(",");
	            		  
	            		  if (key.equals(gtAltKey)) 
	                            altDP= extendedGT.get(key).toString().split(",");
	    			   }
	            	   if (altDP !=null && refDP !=null && refDP.length==2 && altDP.length==2)	  
            		     return (float)Integer.parseInt(altDP[0])/ (Integer.parseInt(refDP[0])+Integer.parseInt(altDP[0]));
	            	   else {
	            		   SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
	            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
	            		  		+ "\t"+"has no strelka AU:CU:GU:TU or correct AU:CU:GU:TU in FORMAT!");
	            		  return -1;  
	            	   }
	    			 }
	    			      			
	    		}

	    	}
	    	else {  //Mutect, Mutect2, Lofreq, Vardict, Varscan, Sniper
	    		String afName="";
//	    		if (callerNameWithHighestPriority.equals("Vardict"))
//	    			System.out.println("found");
	    		if (Integer.bitCount(set)>1) {
		    	   	if (callerNameWithHighestPriority.equals("Mutect")) 
		    	   		afName=callerNameWithHighestPriority+"_FA";
		    	    else 
		    	   		    afName=callerNameWithHighestPriority+"_AF";
		    	   		
		    	    
	    		}
	    		else {
	    			if (callerNameWithHighestPriority.equals("Mutect")) 
		    	   		afName="FA";
		    	    else {
		    	    	if (callerNameWithHighestPriority.equals("Somaticsniper"))
		    	    		afName="DP4";
		    	    	else
		    	    		afName="AF";
		    	   		
		    	    }
	    			
	    		}
	    	   	if (callerNameWithHighestPriority.equals("Lofreq")) {
		    		Map<String,Object> info=variantContext.getAttributes();
		    		for (String key:info.keySet()) {		
		    			if (key.equals(afName))
		    				return Float.parseFloat(info.get(key).toString());	
		    		}
	    	   	}
	    	   	else { //Mutect2 Mutect Vardict VarScan
	    	   		
		    	   		if (variantContext.hasGenotypes()) {
		    	   			boolean foundTumor=false;
		    	   			String tumorSampleName="";
		    	   			Set<String> sampleNames=variantContext.getSampleNames();
		    				for(String sampleName:sampleNames) 
		    					   if (sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG1)||sampleName.toUpperCase().contains(VCFFile.TUMOR_TAG2)) {
		    						  tumorSampleName=sampleName;
		    						  foundTumor=true;  
		    						  break;
		    					   }
		    				if (!foundTumor) {
		    					SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
				            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
				            		  		+ "\t"+"has no Tumor sample in FORMAT for "+ callerNameWithHighestPriority+" calling!");
			    	   			
		    				}
		    				else {
			    			   Genotype gt=variantContext.getGenotype(tumorSampleName);
			    			   Map<String, Object> extendedGT = gt.getExtendedAttributes();
			    			   
			    			   for (String key:extendedGT.keySet()) 
			            		  if (key.equals(afName)) 
			            			//  if (extendedGT.get(afName).toString().toLowerCase().equals("nan"))
		                                 return Float.parseFloat(extendedGT.get(afName).toString());
		                             
			            		  
			    			   SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
				            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
				            		  		+ "\t"+"has no AF or DP in FORMAT for "+callerNameWithHighestPriority+" calling!");	
			    			   
		     
			    			   }
		    	   		}
			    	    else 
			    	    	SomaticCombiner.logger.log(Level.WARNING,variantContext.getContig()+":"+variantContext.getStart()+"\t"+
			            	         variantContext.getAlleles().get(0).getBaseString()+"\t"+variantContext.getAlleles().get(1).getBaseString()
			            		  		+ "\t"+"has no AF in FORMAT for "+ callerNameWithHighestPriority+" calling!");
	    	   				    	   			    	   		
	    	     	
	    	   	}
	    	}
			return -1;
	    		
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

		public int getPriority() {
			return priority;
		}

		public void setPriority(int priority) {
			this.priority = priority;
		}

		
	 }
