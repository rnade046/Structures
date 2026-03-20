package evolutionConservation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class AssessEvolutionConservation {

	public static void main(String[] args) {
		
		String featureBitsFilePrefix = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/evolution/fBits/conservation_motif";
		String evolutionSignificance = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/evolution/fBits/summaryConservation_july2023.tsv";

		int[] modules = new int[] { 82, 85, 94, 239};
		List<String> motifSignificance = new ArrayList<>();

		for(int mod : modules) {
			System.out.println("motif: " + mod);
			FeatureBits fBits = new FeatureBits(featureBitsFilePrefix + mod + ".txt");
			motifSignificance.add(String.valueOf(mod) + "\t" + fBits.getConservationProperties());
		}
		printMotifSignificance(motifSignificance, evolutionSignificance);
	}

	private static void printMotifSignificance(List<String> motifSignificance, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("motif\tConservation\tTrials\tProbSuccess\tFold-change\tp-val\n");
			for(int i=0; i<motifSignificance.size(); i++) {
				out.write(motifSignificance.get(i) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
