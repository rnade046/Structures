package sequenceLength;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

public class Protein {
	
	private String name;
	private HashSet<String> refSeqIds;
	private HashMap<String, Integer> refSeqIdsToLengthMap;
	private int minSeqLength;
	private int maxSeqLength;
	
	public Protein(String n, String[] ids) {
		this.name = n;
		this.refSeqIds = new HashSet<>(Arrays.asList(ids));
		this.refSeqIdsToLengthMap = new HashMap<>();
		this.minSeqLength = Integer.MAX_VALUE;
		this.maxSeqLength = Integer.MIN_VALUE;
	}
	
	public boolean containsId(String id) {
		boolean contains = false;
		if(this.refSeqIds.contains(id)){
			contains = true;
		}
		return contains;
	}
	
	public void setSeqLength(String id, int seqLength) {
		this.refSeqIdsToLengthMap.put(id, seqLength);
		
		if(seqLength < this.minSeqLength) {
			this.minSeqLength = seqLength;
		}
		
		if(seqLength > this.maxSeqLength) {
			this.maxSeqLength = seqLength;
		}
	}
	
	public String getProteinName() {
		return this.name;
	}
	
	public String getIdInfo() {
		String info = "";
		
		for(Entry<String, Integer> entry : refSeqIdsToLengthMap.entrySet()) {
			info += entry.getKey() + "_" + entry.getValue() + "|";
		}
		return info;
	}
	
	public int[] getLengthLimits() {
		return new int[] {minSeqLength, maxSeqLength};
	}
	
	public int getNumberOfIDs() {
		return this.refSeqIds.size();
	}
}