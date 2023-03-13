package structures;

public class End extends SingleStrand{

	String structureType; 
	
	public End(String s) {
		super(s);
		
		structureType = determineStructureType();
	}

	
	private String determineStructureType() {
		String structure = "E" + this.length;
		
		return structure;
	}
	
	public String getStructureType() {
		return structureType;
	}
}
