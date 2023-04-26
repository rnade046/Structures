package checkSequenceStructure;

import java.util.HashMap;
import java.util.HashSet;

public class Occurrence {

	private String protein;
	private String refSeqId;
	private String sequence;
	private char dominantStructure;
	private int structureTypeCount;
	private String structure; 
	private String position;

	public Occurrence(String p, String id, String seq, String struct, int[] pos) {
		this.protein = p;
		this.refSeqId = id;
		this.sequence = seq;
		this.structure = struct; 
		this.dominantStructure = determineDominantStructure();
		this.structureTypeCount = countStructureTypes();
		this.position = String.valueOf(pos[0]) + ".." + String.valueOf(pos[1]);
	}

	private char determineDominantStructure() {

		char s = '*'; 

		int maxOccurrence = 0;

		/* determine occurrence of each character */
		HashMap<Character, Integer> structureCount = new HashMap<>();

		for(char c : this.structure.toCharArray()) {

			if(structureCount.containsKey(c)) {

				int updatedCount = structureCount.get(c) + 1;
				structureCount.put(c, updatedCount);

				if(updatedCount > maxOccurrence) {
					s = c; 
					maxOccurrence = updatedCount;
				}

			} else { 
				structureCount.put(c, 1);
			}
		}
		return s;
	}
	
	private int countStructureTypes() {
		
		HashSet<Character> characterSet = new HashSet<>();
		for(char c : this.structure.toCharArray()) {
			characterSet.add(c);
		}
		return characterSet.size();
	}

	public String getInfo() {
		return new String(this.protein +"\t" +  this.refSeqId + "\t" + this.sequence + "\t" 
				+ this.structure + "\t"+ this.dominantStructure + "\t" + String.valueOf(this.structureTypeCount)+ "\t"+  this.position);
	}
}
