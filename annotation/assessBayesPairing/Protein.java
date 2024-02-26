package assessBayesPairing;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class Protein {

	private String name;
	private HashMap<String, String> jsonFileIdxMap; // fileID = refSeqId;
	private HashSet<String> missedIDs;
	private List<Module> moduleList;
	private List<Double> allScores;
	private HashMap<Integer, Double> moduleSummaryMap; // id = max score obtained

	public Protein(String n) {
		this.name = n;
		this.jsonFileIdxMap = new HashMap<>();
		this.moduleList = new ArrayList<>();
		this.allScores = new ArrayList<>();
		this.missedIDs = new HashSet<>();
	}
	
	public String getProteinName() {
		return this.name;
	}
	
	public HashMap<String, String> getFileIdMap(){
		return this.jsonFileIdxMap;
	}
	
	public List<Module> getModuleList(){
		return this.moduleList;
	}
	
	public HashSet<String> getMissedIDs(){
		return this.missedIDs;
	}

	public void updateJSONmapping(String id, String fileName) {
		this.jsonFileIdxMap.put(fileName, id);
	}
	
	public void addModule(Module m) {
		this.moduleList.add(m);
	}
	
	public void updateMultipleModules(List<Module> m, double threshold) {
		
		for(Module mod: m) {
			/* add module that passes threshold to list */
			if(mod.score > threshold) {
				this.moduleList.add(mod);
			} 
			
			/* add all score to list for percentile calculation */
			this.allScores.add(mod.getScore());
		}
	}
	
	public void addMissedId(String id) {
		this.missedIDs.add(id);
	}

	public void summarizeModules() {
		
		HashMap<Integer, Double> mMap = new HashMap<>();
		
		for(Module m : moduleList) {
			
			int id = m.getID();
			double score = m.getScore();
			
			if(mMap.containsKey(id)) {
				double currentMaxScore = mMap.get(id);
				if(score > currentMaxScore) {
					mMap.put(id, score);
				}
			} else {
				mMap.put(id, score);
			}
		}
		this.moduleSummaryMap = mMap;
	}
	
	public HashMap<Integer, Double> getModuleSummaryMap(){
		return this.moduleSummaryMap;
	}
	
	public List<Double> getAllScores(){
		return this.allScores;
	}
}
