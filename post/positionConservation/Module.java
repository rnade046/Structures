package positionConservation;

public class Module {

	int seqLength;
	double score;
	int[] positions;
	
	public Module(int sLength, double s, String[] p) {
		seqLength = sLength;
		score = s;
		positions = setPositions(p);
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
	
	private int[] setPositions(String[] p) {
		
		int[] pos3 = new int[p.length];
		
		for(int j=0; j<p.length; j++) {
			pos3[j] = Integer.parseInt(p[j]);
		}
		
		return pos3;
	}

}
