package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;

import org.json.JSONArray;
import org.json.JSONObject;

public class assessPositionConservationInBins {

	private static int bins = 50;

	public static void main(String[] args) throws FileNotFoundException, IOException {

		//Properties params = new Properties();
		//params.load(new FileInputStream(args[0]));	

		//String wd = params.getProperty("working_directory");
		String wd = args[0];
		String jsonPath = args[1];
		String significantModulesFile = wd + args[2];

		Files.createDirectories(Paths.get(wd + "positionCons/"));
		String positionConservationFilePrefix = wd + "positionCons/positionConservation_m";

		/* search each module for position conservation */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(significantModulesFile))));		

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				/* obtain module info */
				String[] info = line.split("\t");
				String module = info[0];

				if(info[1].equals("HL")) {
					System.out.println("testing module: " + module + " | #prots : " + info[4]);

					/* initialize motif positions list and considered sequences for normalization */
					int[] motifPositions = new int[bins];
					int[] consideredSequence = new int[bins];
					double[] normalizedPositions = new double[bins];

					for(String protein : info[5].split("\\|")) {

						System.out.println(protein.split("\\_",2)[0] + " : ");

						String[] ids = protein.split("\\_",2)[1].substring(1, protein.split("\\_", 2)[1].length()-1).split(",");

						Double score = 0.0;
						int seqLength = 0;
						int[] positions = null;
						String finalIdentifiers = "";

						for(String i: ids) { // refSeqId=jsonPath

							String[] identifiers = i.split("=");

							/* search for module with highest score within each sequence corresponding to refSeqId */
							if(identifiers.length > 1) {
								File f = new File(jsonPath + identifiers[1]);
								if(f.exists()) {

									/* obtain information for current refSeqId*/
									Module mod = searchJsonForMaxModuleScore(jsonPath + identifiers[1], identifiers[0], module, score);

									/* update information if values array is not empty (i.e. corresponds to modules with > score) */
									if(mod != null) {
										score = mod.getScore();
										positions = mod.getPositions();
										seqLength = mod.getSeqLength() - 100;  // Subtract 100 CDS
										finalIdentifiers = identifiers[0];
									}
								} 
							}
						} 

						/* determine position of current module */
						if(positions != null) {
							System.out.println(finalIdentifiers + " | score: " + score);

//							if(Math.floor(seqLength / bins) >= motifPositions.length) {
								int[][] binnedPositions = formatSequenceInBins(seqLength);
								motifPositions = increasePositionCount(positions, binnedPositions, motifPositions);
								consideredSequence = updateConsideredSequencesForNormalization(positions, binnedPositions, consideredSequence);
//							} 
						}
					}
					/* normalize positions */
					normalizedPositions = normalizePositionsByConsideredPositions(motifPositions, consideredSequence, normalizedPositions);
					printNormalizedMotifPosition(normalizedPositions, positionConservationFilePrefix + module + ".tsv");
				}
				System.out.println();
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static Module searchJsonForMaxModuleScore(String inputFile, String id, String module, Double maxScore) {

		Module mod = null;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();

			JSONObject rootObject = new JSONObject(line); // Parse the JSON to a JSONObject
			int seqLength = rootObject.getJSONArray("input").get(0).toString().length();

			if(rootObject.getJSONObject("all_hits").has(id)){

				JSONObject a = rootObject.getJSONObject("all_hits").getJSONObject(id); // {270}
				//				System.out.print(id + ": ");

				if(a.has(module)) {

					JSONArray m = a.getJSONArray(module); 
					for(int i=0; i < m.length(); i++) {

						/* each instance of a module */
						JSONArray elements = m.getJSONArray(i);
						//String seq = elements.get(0).toString();
						String[] pos = elements.get(1).toString().split("[,\\[\\]]");
						String[] pos2 = Arrays.copyOfRange(pos, 2, pos.length);
						double score = Double.parseDouble(elements.get(2).toString());

						if(score > 4) { // positive score
							if(!pos2[0].isEmpty()) {
								if(Integer.parseInt(pos2[0]) > 100) { // ignore models in the CDS

									if(score > maxScore) {
										mod = new Module(seqLength, score, pos2);
									}	
								}
							}
						}
					}
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return mod;
	}

	private static int[][] formatSequenceInBins(int seqLength) {

		int[][] sequenceBins = new int[bins][2];

		int[] binSizes = new int[bins];
		int[] binSizesReodered = new int[bins];

		int minimumBinSize = (int) Math.floor(seqLength / (double) bins); 
		int binRemainder = seqLength % bins; 

		/* determine size of bins (adding remainder to minimum bin sizes) */
		for(int i=0; i<bins; i++) {

			// determining indexes for unequal bins : https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
			int binStart = i*minimumBinSize + Math.min(i, binRemainder);
			int binEnd = (i+1) * minimumBinSize + Math.min(i+1, binRemainder);
			int binRange = binEnd - binStart;

			binSizes[i] = binRange;
		}

		/* re-ordering bin sizes to pad edges */
		int countEven = 0;
		int countOdd = bins-1;
		for(int i=0; i<bins; i++) {
			if(i%2 ==0) {
				binSizesReodered[countEven] = binSizes[i];
				countEven++;
			} else { 
				binSizesReodered[countOdd] = binSizes[i];
				countOdd--;
			}
		}

		/* bin positions (0 = start position, 1 = end position) */
		int currentIdx = 100;
		for(int i=0; i<bins; i++) {

			sequenceBins[i][0] = currentIdx; // start
			sequenceBins[i][1] = currentIdx + binSizesReodered[i]; // end
			currentIdx += binSizesReodered[i];
		}
		return sequenceBins;
	}

	private static int[] increasePositionCount(int[] positions, int[][] binnedPositions, int[] motifPositions) { 

		boolean motifFound = false;

		/* iterate over each bin */
		for(int bin=0; bin<binnedPositions.length; bin++) {

			int start = binnedPositions[bin][0];
			int end = binnedPositions[bin][1];

			if(positions[0] >= start && positions[positions.length-1] <= end) {
				motifPositions[bin] += 1;
				motifFound = true;
				break;  // stop search
			} 
		}

		if(!motifFound) { 
			System.out.println("module accross multiple bins");
		}
		return motifPositions;
	}

	private static int[] updateConsideredSequencesForNormalization(int[] positions, int[][] binnedPositions, int[] consideredSequences){

		/* iterate over each bin */
		for(int bin=0; bin<binnedPositions.length; bin++) {

			/* determine number of motifs tested for each subsequence (i.e. bins)  */
			int currentLength = (binnedPositions[bin][1] - binnedPositions[bin][0]) - positions.length + 1;
			consideredSequences[bin] += currentLength;
		}
		return consideredSequences;
	}

	private static double[] normalizePositionsByConsideredPositions(int[] motifPositions, int[] consideredSeq, double[] normalizedPositions) {

		for(int bin=0; bin<normalizedPositions.length; bin++) {
			normalizedPositions[bin] = motifPositions[bin] / (double) consideredSeq[bin];
		}
		return normalizedPositions;
	}

	private static void printNormalizedMotifPosition(double[] normalizedPositions, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<normalizedPositions.length; i++) {
				out.write((i+1) + "\t" + normalizedPositions[i] + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
