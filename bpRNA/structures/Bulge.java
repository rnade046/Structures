package structures;

public class Bulge extends SingleStrand{

	String structureType; 
	
	public Bulge(String s) {
		super(s);
		
		structureType = determineStructureType();
	}

	
	private String determineStructureType() {
		String structure = "B" + this.length;
		
		return structure;
	}
	
	public String getStructureType() {
		return structureType;
	}
}
