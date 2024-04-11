package positionConservation;

public class Module {

	private int seqLength;
	private double score;
	private int[] positions;
	private String identifiers;


	public Module(int sLength, double s) {
		this.seqLength = sLength;
		this.score = s;
		this.positions = new int[2];
	}	

	public Module(int sLength, double s, String[] p, boolean utr) {
		this.seqLength = sLength;
		this.score = s;

		if(utr) {
			this.positions = setUTRPositions(p);
		} else {
			this.positions = setPositions(p);
		}

	}

	public int getSeqLength() {
		return this.seqLength;
	}

	public double getScore() {
		return this.score;
	}

	public int[] getPositions(){
		return this.positions;
	}

	public void setIdentifiers(String id) {
		this.identifiers = id;
	}

	public String getIdentifiers() {
		return this.identifiers;
	}

	private int[] setPositions(String[] p) {

		int[] pos3 = new int[p.length];

		for(int j=0; j<p.length; j++) {
			if(!p[j].isEmpty()) {
				pos3[j] = Integer.parseInt(p[j]);
			}
		}

		return pos3;
	}

	private int[] setUTRPositions(String[] p) {

		int[] pos3 = new int[p.length];

		for(int j=0; j<p.length; j++) {
			if(!p[j].isEmpty()) {
				pos3[j] = Integer.parseInt(p[j]) -100;
			}
		}

		return pos3;
	}

}
