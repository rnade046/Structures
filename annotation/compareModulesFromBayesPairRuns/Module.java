package compareModulesFromBayesPairRuns;

import java.util.HashSet;

public class Module {

	private String moduleId;

	private HashSet<String> newProteins;
	private HashSet<String> lostProteins;
	private HashSet<String> keptProteins;

	public Module(String id) {
		this.moduleId = id;

		this.newProteins = new HashSet<>();
		this.lostProteins = new HashSet<>();
		this.keptProteins = new HashSet<>();
	}

	public String getModuleId() {
		return this.moduleId;
	}

	public HashSet<String> getNewProteins(){
		return this.newProteins;
	}

	public HashSet<String> getLostProteins(){
		return this.lostProteins;
	}

	public HashSet<String> getKeptProteins(){
		return this.keptProteins;
	}
}
