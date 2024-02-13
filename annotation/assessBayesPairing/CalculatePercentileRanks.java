package assessBayesPairing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

public class CalculatePercentileRanks {

	public static void determinePercentiles(List<Protein> proteinList, String outputFile) {
		// test
		System.out.println("this is a test");
		/* obtain scores in ascending order */
		List<Double> scores = obtainScores(proteinList);

		/* assign percentiles to scores */ 
		HashMap<Double, List<Double>> percentiles = evaluateScoreForPercentiles(scores);

		/* print output */
		printPercentiles(percentiles, outputFile);
	}

	/**
	 * Obtain all scores in ascending order from all proteins (sequences searched with BayesPairing)
	 * @param proteinList	List<Protein> - Protein objects contains Module list
	 * @return scores 		List<Double> - scores in ascending order 
	 */
	public static List<Double> obtainScores(List<Protein> proteinList) {

		ArrayList<Double> scores = new ArrayList<>();

		/* search all proteins */
		for(Protein prot: proteinList) {

			/* obtain scores for all modules */
			scores.addAll(prot.getAllScores());
		}

		/* sort scores in ascending order */
		Collections.sort(scores);

		return scores;
	}
	/**
	 * Assign percentile to each score
	 * 
	 * @param scores		List<Double> - list of all scores (from modules identified w/ BayesParing)
	 * @return percentiles	HashMap<Double, List<Double>> - map of {percentile (#.##) = list(scores)}  	
	 */
	public static HashMap<Double, List<Double>> evaluateScoreForPercentiles(List<Double> scores) {

		HashMap<Double, List<Double>> percentiles = new HashMap<>();

		for(int i=1; i<scores.size(); i++) {

			//  ceiling; to prevent 0 probabilities
			double p = Math.ceil(((i-1) / (double) scores.size()) * 100.00);
			//			System.out.println(i + "|" + p + "|" + ((i-1)/ (double) scores.size()));
			if(percentiles.containsKey(p)) { 
				percentiles.get(p).add(scores.get(i));
			} else { 
				ArrayList<Double> l = new ArrayList<>();
				l.add(scores.get(i));
				percentiles.put(p, l);
			}
		}
		return percentiles;
	}

	/**
	 * Print percentiles
	 * Percentile (#.##) | MinScore_MaxScore | List of all scores |
	 * 
	 * @param percentiles
	 * @param outputFile
	 */
	public static void printPercentiles(HashMap<Double, List<Double>> percentiles, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(Entry<Double, List<Double>> entry : percentiles.entrySet()) {

				List<Double> scores = entry.getValue();
				Collections.sort(scores);

				out.write(entry.getKey() + "\t" + scores.get(0) + "_" + scores.get(scores.size()-1) + "\t" + scores + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
