package positionConservation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

public class Segment {

	private int[] motifPositions;
	private int[] consideredSequence;
	private double[] normalizedPositions;

	public Segment(int length) {
		this.motifPositions = new int[length];
		this.consideredSequence = new int[length];
		this.normalizedPositions = new double[length];
	}

	public void utrIncreasePositionCount(Module mod) { 

		/* find corresponding position for each nucleotide position */
		int halfLengthCeil = (int) Math.ceil(mod.getSeqLength() / (double) 2);
		int halfLengthFloor = (int) Math.floor(mod.getSeqLength() / (double) 2);

		int[] positions = mod.getPositions();
		System.out.println("pos: " + Arrays.toString(mod.getPositions()) + " | hlC= " + halfLengthCeil + " | hlF= " + halfLengthFloor);

		for(int p=0; p<positions.length; p++) {

			int currentPosition = positions[p];

			if((p==0 && currentPosition==0) || currentPosition != 0) {
				if(currentPosition < halfLengthCeil) {
					this.motifPositions[currentPosition]++;
				} else {
					this.motifPositions[motifPositions.length - (mod.getSeqLength() - currentPosition)]++;
				}
			}
		}

		/* adjust considered sequences */
		for(int i=0; i<= halfLengthCeil; i++) {
			this.consideredSequence[i]++;
		}
		for(int i=consideredSequence.length-1; i>=halfLengthFloor; i--) {
			this.consideredSequence[i]++;
		}
	}

	public void cdsIncreasePositionCount(Module mod) {

		/* increase position count */

		System.out.println("pos: " + Arrays.toString(mod.getPositions()) + " | l= " + mod.getSeqLength());
		int[] positions = mod.getPositions();
	
		for(int p=0; p<positions.length; p++) {
			if((p==0 && positions[p]== 0) || positions[p] != 0) {
				this.motifPositions[positions[p]]++;
			}
		}

		/* adjust considered sequences (note : always 100 for CDS) */
		for(int i=0; i< mod.getSeqLength(); i++) {
			this.consideredSequence[i]++;
		}
	}

	public void normalizePositionsByConsideredPositions() {

		for(int i=0; i<this.normalizedPositions.length; i++) {
			this.normalizedPositions[i] = this.motifPositions[i] / (double) this.consideredSequence[i];
		}
	}

	public void printMotifPosition(String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<this.normalizedPositions.length; i++) {
				out.write((i+1) + "\t" + this.motifPositions[i] + "\t" + this.consideredSequence[i] + "\t" + this.normalizedPositions[i] + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
