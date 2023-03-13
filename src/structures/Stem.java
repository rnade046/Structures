package structures;

public class Stem {

	String id;
	int[] stem1Bounds;
	int[] stem2Bounds;
	String stem1Seq;
	String stem2Seq; 
	int length;
	String structureType;
	
	public Stem(String s) {
		String[] elements = s.split("\\s+");
		
		id = elements[0];
		
		stem1Bounds = formatBounds(elements[1]);
		stem2Bounds = formatBounds(elements[3]);
		
		stem1Seq = elements[2].split("\"")[1];
		stem2Seq = elements[4].split("\"")[1];
		
		length = stem1Bounds[1] - stem1Bounds[0] + 1;
		
		structureType = determineStructureType();
	}
	
	private static int[] formatBounds(String elements) {
	
		String[] bounds = elements.split("\\..");
		int[] boundsFormated = new int[bounds.length];
		
		for(int i=0; i<bounds.length; i++) {
			boundsFormated[i] = Integer.parseInt(bounds[i]);
		}
		
		return boundsFormated;
	}
	
	private String determineStructureType() {
		return new String("S" + length);
	}
	
	public String getStructureType() {
		return structureType;
	}
	
	
	
 }
