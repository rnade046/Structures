package info;

import java.util.HashSet;
import java.util.Set;

public class Module {

	private String name;
	private int numberProteins;
	private double pval;
	private HashSet<String> proteinSet;
	private String type;
	private String atlasID;
	
	public Module(String n, int nProts, double p) {
		this.name = n;
		this.numberProteins = nProts;
		this.pval = p;
		this.proteinSet = new HashSet<>();
	}
	
	public Set<String> getProteins() {
		return this.proteinSet;
	}
	
	public void setMetaData(String aID, String t) {
		this.atlasID = aID;
		this.type = t;
	}
	
	public void setProteins(String[] proteins) {
		for(String prot: proteins) {
			this.proteinSet.add(prot.split("\\_")[0]);
		}
	}
	
	public boolean containsProtein(String protein) {
		boolean contains = false;
		if(this.proteinSet.contains(protein)) {
			contains = true;
		}
		return contains;
	}
	
	public String[] getModuleSummary() {
		return new String[] {this.name, this.type, this.atlasID,
				String.valueOf(this.pval), String.valueOf(this.numberProteins)};
	}
}
