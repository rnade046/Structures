package assessBayesPairing;

public class Module {

	int id;
	String sequence; 
	double score;
	String[] positions;
	
	public Module(int i, String seq, double s, String[] p) {
		id = i;
		sequence = seq;
		score = s;
		positions = p;
	}
	
	public int getID() {
		return this.id;
	}
	
	public String getSequence() {
		return this.sequence;
	}
	
	public double getScore() {
		return this.score;
	}
	
	public String[] getPositions(){
		return this.positions;
	}

}
