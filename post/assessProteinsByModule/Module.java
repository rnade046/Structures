package assessProteinsByModule;

import java.util.HashSet;

public class Module {

	private String name;
	private int numberProteins;
	private HashSet<String> proteinSet;
	
	public Module(String n, int nProts) {
		this.name = n;
		this.numberProteins = nProts;
		this.proteinSet = new HashSet<>();
	}
	
	public String getProteinName() {
		return this.name;
	}
	
	public int getNumberOfProteins() {
		return this.numberProteins;
	}
	
	public void setProteins(String[] proteins) {
		for(String prot: proteins) {
			this.proteinSet.add(prot.split("\\_")[0]);
		}
	}
	
	public HashSet<String> getProteins(){
		return this.proteinSet;
	}
	
	public boolean containsProtein(String protein) {
		boolean contains = false;
		if(this.proteinSet.contains(protein)) {
			contains = true;
		}
		return contains;
	}
}
