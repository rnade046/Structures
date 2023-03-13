package structures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class MultiLoop {

	String id;
	List<int[]> bounds;
	List<String> seqs;
	List<Integer> lengths;
	String structureType;
	
	public MultiLoop(List<String> s) {
		
		bounds = new ArrayList<>();
		seqs = new ArrayList<>();
		lengths = new ArrayList<>();
		
		id = s.get(0).split("\\.")[0];
		
		for(int i=0; i<s.size(); i++) {
			
			String[] elements = s.get(i).split("\\s+");	
			int[] b = formatBounds(elements[1]);
			bounds.add(b);
			
			
			String[] seq = elements[2].split("\"");
			if(seq.length > 0) {
				seqs.add(seq[0]);
			} else {
				seqs.add("");
			}
			
			lengths.add(b[1] - b[0] + 1);
			
		}
		
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
		String structure = "M";
		
		Collections.sort(lengths);
		for(int i=0; i<lengths.size(); i++) {
			if(i<lengths.size()-1) {
				structure += lengths.get(i) + ".";
			} else {
				structure += lengths.get(i);
			}
		}
		
		return structure;
	}
	
	public String getStructureType() {
		return structureType;
	}
 }
