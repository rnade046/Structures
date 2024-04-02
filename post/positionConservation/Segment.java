package positionConservation;

public class Segment {

	int[] motifPositions;
	int[] consideredSequence;
	double[] normalizedPositions;
	
	public Segment(int length) {
		this.motifPositions = new int[length];
		this.consideredSequence = new int[length];
		this.normalizedPositions = new double[length];
	}
	
}
