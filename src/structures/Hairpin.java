package structures;

public class Hairpin extends SingleStrand{

	String structureType; 
	
	public Hairpin(String s) {
		super(s);
		
		structureType = determineStructureType();
	}

	
	private String determineStructureType() {
		String structure = "H" + this.length;
		
		return structure;
	}
	
	public String getStructureType() {
		return structureType;
	}
}

