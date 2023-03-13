package structures;

public class SingleStrand {

	String id;
	int[] bounds;
	String seq;
	int length;
	
	public SingleStrand(String s) {
		String[] elements = s.split("\\s+");
		
		id = elements[0];
		
		bounds = formatBounds(elements[1]);
		
		seq = elements[2].split("\"")[0];
		
		length = bounds[1] - bounds[0] + 1;
	}
	
	private static int[] formatBounds(String elements) {
		
		String[] bounds = elements.split("\\..");
		int[] boundsFormated = new int[bounds.length];
		
		for(int i=0; i<bounds.length; i++) {
			boundsFormated[i] = Integer.parseInt(bounds[i]);
		}
		
		return boundsFormated;
	}
	
}
