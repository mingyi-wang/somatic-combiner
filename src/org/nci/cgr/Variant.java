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
	protected VariantContext variantContext;
	    // private String variantKey;
	    protected String caller;
	    protected Byte set;
	    protected int priority;
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

		
	 }
