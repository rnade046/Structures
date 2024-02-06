package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

public class ScorePercentiles {

	public static void assignScorePercentile(String annotationFile, String percentileFile, String updatedAnnotationFile, double rescaleThreshold) {

		/* load percentiles from file; where i = percentile, j(1) >= score, j(2) <= score */
		List<Double> scores = loadScoresByPercentileThreshold(percentileFile, rescaleThreshold);

		/* scale scores into percentiles */
		HashMap<Double, HashSet<Double>> percentiles = evaluateScoreForPercentiles(scores);
		HashMap<Double, Double[]> percentileBounds = obtainPercentileBounds(percentiles);
		/* print rescaled percentiles */
		printPercentiles(percentiles, percentileFile, rescaleThreshold);

		/* update scores within annotation file */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(updatedAnnotationFile)));

			String line = in.readLine();
			out.write(line + "\n"); //header

			line = in.readLine();

			while(line!=null) {

				String[] col = line.split("\t"); // [0] = module, [1] = #proteins, [2]=list proteins(prot1_score1|..|)
				if(col.length > 2) {
					out.write(col[0] + "\t" + col[1] + "\t");

					String[] prots=col[2].split("\\|");
					for(String p : prots) {

						String[] info = p.split("\\_"); // [0] = gene_name, [1] = score

						/* determine percentile */
						double percentile = determineScorePercentile(Double.parseDouble(info[1]), percentileBounds);
						System.out.println(info[0] + "_" + info[1] + " _" + percentile);
						out.write(info[0] + "_" + percentile + "|");
					}
					out.write("\n");
					out.flush();
				}
				line = in.readLine();
			}
			in.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static List<Double> loadScoresByPercentileThreshold(String inputFile, double threshold){

		/* i = percentile, j(1) >= score, j(2) <= score */
		List<Double> scores = new ArrayList<>();

		/* load scores from percentiles file */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();
			while(line!=null) {

				String[] col = line.split("\t"); //[0] = percentile, [1] = min_max, [2] = {score1, score2, ..., scoreN}

				int p = (int) Double.parseDouble(col[0]);

				if(p >=threshold*100) {
					String[] values = col[2].substring(1, col[2].length()-1).split("[\\s,]+");

					for(String v: values) {
						scores.add(Double.valueOf(v));
					}
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return scores;
	}

	/**
	 * Assign percentile to each score
	 * 
	 * @param scores		List<Double> - list of all scores (from modules identified w/ BayesParing)
	 * @return percentiles	HashMap<Double, List<Double>> - map of {percentile (#.##) = list(scores)}  	
	 */
	public static HashMap<Double, HashSet<Double>> evaluateScoreForPercentiles(List<Double> scores) {

		HashMap<Double, HashSet<Double>> percentiles = new HashMap<>();

		for(int i=1; i<scores.size(); i++) {

			//  ceiling; to prevent 0 probabilities
			double p = Math.ceil(((i-1) / (double) scores.size()) * 100.00);
			//			System.out.println(i + "|" + p + "|" + ((i-1)/ (double) scores.size()));
			if(percentiles.containsKey(p)) { 
				percentiles.get(p).add(scores.get(i));
			} else { 
				HashSet<Double> l = new HashSet<>();
				l.add(scores.get(i));
				percentiles.put(p, l);
			}
		}
		return percentiles;
	}

	public static double determineScorePercentile(double score, HashMap<Double, Double[]> percentiles)   {

		double percentile = 0.0;

		for(Entry<Double, Double[]> e: percentiles.entrySet()) {

			Double[] bounds = e.getValue();
			
			if(score >= bounds[0] && score <= bounds[1]) {
				percentile = e.getKey();
				break;
			}
		}
		return percentile;
	}

	public static HashMap<Double, Double[]> obtainPercentileBounds(HashMap<Double, HashSet<Double>> percentiles) {

		HashMap<Double, Double[]> percentileBounds = new HashMap<>();
		for(Entry<Double, HashSet<Double>> entry : percentiles.entrySet()) {

			List<Double> scores = new ArrayList<Double>(entry.getValue());
			Collections.sort(scores);

			percentileBounds.put(entry.getKey(), new Double[] {scores.get(0), scores.get(scores.size()-1)});
		}
		return percentileBounds;
	}

	/**
	 * Print percentiles
	 * Percentile (#.##) | MinScore_MaxScore | List of all scores |
	 * 
	 * @param percentiles
	 * @param outputFile
	 */
	public static void printPercentiles(HashMap<Double, HashSet<Double>> percentiles, String outputFile, double threshold) {

		try {
			String file = outputFile.substring(0, outputFile.length()-4);
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(file + "_rescaled_"+ threshold +".tsv")));

			for(Entry<Double, HashSet<Double>> entry : percentiles.entrySet()) {

				List<Double> scores = new ArrayList<Double>(entry.getValue());
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
