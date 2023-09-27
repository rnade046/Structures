package assessBayesPairing;

public class JSON {

	private String jsonFile;
	private String refSeqId;
	private String protein;
	private int seqLength;
	private boolean exists;
	
	public JSON(String json, String id, int l) {
		this.jsonFile = json;
		this.refSeqId = id;
		this.seqLength = l; 
		this.exists = false;
	}
	
	public void fileExists() {
		this.exists = true;
	}
	
	public void setProtein(String p) {
		this.protein = p;
	}
	
	public String[] getEntry() {
		return new String[] {this.jsonFile, this.refSeqId, this.protein, String.valueOf(this.seqLength), String.valueOf(this.exists)};
	}
}
