package input;

public class Protein {
	
	private String proteinName;
	private int safeDomain;
	private int nmfRegions;
	
	
	public Protein(String name, int domain, int region) {
		this.proteinName = name;
		this.safeDomain = domain;
		this.nmfRegions = region;
	}
	
	public String getName() {
		return this.proteinName;
	}
	
	public int getSafeDomain() {
		return this.safeDomain;
	}
	
	public int getNMFregion() {
		return this.nmfRegions;
	}
}
