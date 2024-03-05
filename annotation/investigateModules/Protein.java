package investigateModules;

public class Protein {

	private String name;
	private int matrixIndex;
	private int degrees;
	
	public Protein(String n, int idx) {
		this.name = n;
		this.matrixIndex = idx;
		this.degrees = 0;
	}
	
	public String getName() {
		return this.name;
	}
	
	public void setDegrees(int d) {
		this.degrees = d;
	}
	
	public void increaseDegrees() {
		this.degrees++;
	}
	
	public int getIdx() {
		return this.matrixIndex;
	}
	
	public int getDegree() {
		return this.degrees;
	}
}
