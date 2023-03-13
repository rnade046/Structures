package structures;

public class InnerLoop {

	String id;
	int[] loop1Bounds;
	int[] loop2Bounds;
	String loop1Seq;
	String loop2Seq; 
	int length1;
	int length2;
	String structureType;
	
	public InnerLoop(String s1, String s2) {
		String[] elements1 = s1.split("\\s+");
		String[] elements2 = s2.split("\\s+");
		
		id = elements1[0].split("\\.")[0];
		
		loop1Bounds = formatBounds(elements1[1]);
		loop2Bounds = formatBounds(elements2[1]);
		
		loop1Seq = elements1[2].split("\"")[0];
		loop2Seq = elements2[4].split("\"")[0];
		
		length1 = loop1Bounds[1] - loop1Bounds[0] + 1;
		length2 = loop2Bounds[1] - loop2Bounds[0] + 1;
		
		structureType = determineStructureType();
	}
	
	private int[] formatBounds(String elements) {
	
		String[] bounds = elements.split("\\..");
		int[] boundsFormated = new int[bounds.length];
		
		for(int i=0; i<bounds.length; i++) {
			boundsFormated[i] = Integer.parseInt(bounds[i]);
		}
		
		return boundsFormated;
	}
	
	private String determineStructureType() {
		String structure = "I";
		
		if(this.length1 < this.length2) {
			structure += length1 + "." + length2;
			
		} else { 
			structure += length2 + "." + length1;
		}
		
		return structure;
	}
	
	public String getStructureType() {
		return structureType;
	}
 }
