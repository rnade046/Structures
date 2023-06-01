package input;

import java.util.HashMap;
import java.util.HashSet;

public class MotifStructures {

	// STRUCTURE = List<PROTEIN> 
	private HashMap<Character, HashSet<String>> structures; 
	
	public MotifStructures() {
		
		/* initialize structures */
		HashMap<Character, HashSet<String>> struct = new HashMap<>();
		char[] dominantStructures = new char[] {'B', 'E', 'H', 'I', 'M', 'S', 'X'}; 
		for(Character s: dominantStructures) {
			struct.put(s, new HashSet<String>());
		}
		
		this.structures = struct;
	}
	
	public void updateStructures(Character s, String protein) {
		
		HashSet<String> currentStructures = this.structures.get(s);
		currentStructures.add(protein);
		this.structures.put(s, currentStructures);
	}
	
	public HashMap<Character, HashSet<String>> getStructures(){
		return this.structures;
	}
}
