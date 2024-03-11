package assessMembership;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Module {

	private String name;
	private int numberProteins;
	private List<String> proteinSet;
	private List<Double> proteinScores;
	private HashMap<String, Integer> proteinIdxMap;
	
	public Module(String n, String[] proteinInfo) {
		this.name = n;
		this.numberProteins = proteinInfo.length;
		this.proteinSet = setProteins(proteinInfo);
		this.proteinScores = setScores(proteinInfo);
		this.proteinIdxMap = setProteinIdxMap();
	}
	
	public Module(String n, List<String> prots) {
		this.name = n;
		this.proteinSet = prots;
		this.numberProteins = prots.size();
	}
	
	public String getProteinName() {
		return this.name;
	}
	
	public int getNumberOfProteins() {
		return this.numberProteins;
	}
	
	private List<String> setProteins(String[] proteinInfo) {
		
		List<String> proteins = new ArrayList<>();
		
		for(String prot: proteinInfo) {
			proteins.add(prot.split("\\_")[0]);
		}
		return proteins;
	}
	
	private List<Double> setScores(String[] proteins) {
		
		List<Double> scores = new ArrayList<>();
		
		for(String prot: proteins) {
			scores.add(Double.parseDouble(prot.split("\\_")[1]));
		}
		return scores;
	}
	
	public List<String> getProteins(){
		return this.proteinSet;
	}
	
	private HashMap<String, Integer> setProteinIdxMap(){
		
		HashMap<String, Integer> proteinIdxMap = new HashMap<>();
		
		for(int i=0; i<this.proteinSet.size(); i++) {
			proteinIdxMap.put(this.proteinSet.get(i), i);
		}
		
		return proteinIdxMap;
	}
	
	public List<Double> getProteinScores(){
		return this.proteinScores;
	}
	
	public HashMap<String, Integer> getProteinIdxMap(){
		return this.proteinIdxMap;
	}
	
	public boolean containsProtein(String protein) {
		boolean contains = false;
		if(this.proteinSet.contains(protein)) {
			contains = true;
		}
		return contains;
	}
}
