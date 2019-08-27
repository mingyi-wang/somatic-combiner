package org.nci.cgr;

public class Caller {
	private String name;
	private String description;
	private String symbol;
	private Byte set;
    private int priority;
    private String filePath;
    private String type;
    private VCFFile vcfFile;
    private String callerName;
	public Caller(String name, String description, String symbol, Byte set, int priority,String type,String callerName) {
		super();
		this.name = name;
		this.description = description;
		this.symbol = symbol;
		this.set = set;
		this.priority = priority;
		this.filePath=null;
		this.type=type;
		this.vcfFile=null;
		this.callerName=callerName;
	}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getDescription() {
		return description;
	}
	public void setDescription(String description) {
		this.description = description;
	}
	public String getSymbol() {
		return symbol;
	}
	public void setSymbol(String symbol) {
		this.symbol = symbol;
	}
	public Byte getSet() {
		return set;
	}
	public void setSet(Byte set) {
		this.set = set;
	}
	public int getPriority() {
		return priority;
	}
	public void setPriority(int priority) {
		this.priority = priority;
	}

	public String getFilePath() {
		return filePath;
	}

	public void setFilePath(String filePath) {
		this.filePath = filePath;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public VCFFile getVcfFile() {
		return vcfFile;
	}

	public void setVcfFile(VCFFile vcfFile) {
		this.vcfFile = vcfFile;
	}

	public String getCallerName() {
		return callerName;
	}

	public void setCallerName(String callerName) {
		this.callerName = callerName;
	}
    
	

}
