package structures;

public class ExternalLoop extends SingleStrand{

	String structureType; 
	
	public ExternalLoop(String s) {
		super(s);
		
		structureType = determineStructureType();
	}

	
	private String determineStructureType() {
		String structure = "X" + this.length;
		
		return structure;
	}
	
	public String getStructureType() {
		return structureType;
	}
}
