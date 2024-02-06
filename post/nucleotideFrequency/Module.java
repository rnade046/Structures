package nucleotideFrequency;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class Module {

	private HashSet<String> proteins;
	private HashSet<String> refSeqIds;
	private HashSet<String> jsonPaths;

	public Module(String[] annotatedProteins, HashMap<String, String[]> idMap, HashMap<String, String> jsonIdxMap) { 
		
		this.proteins = setProteins(annotatedProteins);
		this.refSeqIds = setRefSeqIds(idMap);
		this.jsonPaths = setJSONpaths(jsonIdxMap);

	}

	private HashSet<String> setProteins(String[] proteins) { 
		HashSet<String> protSet = new HashSet<String>();

		for(String p: proteins) {
			protSet.add(p.split("\\_")[0]);
		}
		return protSet;
	}

	private HashSet<String> setRefSeqIds(HashMap<String, String[]> idMap){
		HashSet<String> idSet = new HashSet<>();
		
		for(String p: this.proteins) {
			if(idMap.containsKey(p)) {
				idSet.addAll(Arrays.asList(idMap.get(p)));
			}
		}
		return idSet;
	}
	
	private HashSet<String> setJSONpaths(HashMap<String, String> jsonIdxMap){
		HashSet<String> pathSet = new HashSet<>();
		
		for(String id: this.refSeqIds) {
			if(jsonIdxMap.containsKey(id)) {
				pathSet.add(jsonIdxMap.get(id));
			}
		}
		return pathSet;
	}
	
	public HashSet<String> getJsonPaths() {
		return this.jsonPaths;
	}
}
