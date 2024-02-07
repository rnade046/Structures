package info;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Protein {
	
	private String name;
	private HashMap<String, String> ids; // refSeqId = jsonIdx
	
	public Protein(String n) {
		this.name = n;
		this.ids = new HashMap<>();
	}

	public String getName() {
		return this.name;
	}
	
	public Set<Map.Entry<String, String>> getIds(){
		return this.ids.entrySet();
	}
	
	public Set<String> getRefSeqIds(){
		return this.ids.keySet();
	}
	
	public Collection<String> getJSONpaths(){
		return this.ids.values();
	}
	
	public void setIds(String[] refSeqIds) {
		for(String i: refSeqIds) {
			this.ids.put(i, "");
		}
	}
	
	public void setJSONPath(String refSeqId ,String jsonFile) {
		this.ids.put(refSeqId, jsonFile);
	}
}
