package sampling;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class ProteinAnnotations {

	private int lowerBoundToSample;
	private int upperBoundToSample;

	private HashSet<String> proteinSet;

	public ProteinAnnotations(int _lowerBoundToSample, int _upperBoundToSample, HashSet<String> _proteinSet) {
		this.lowerBoundToSample = _lowerBoundToSample;
		this.upperBoundToSample = _upperBoundToSample;

		this.proteinSet = _proteinSet;
	}

	/**
	 * Load annotation list and store in HashMap proteins and their number of occurrence in annotation list.
	 * @param annotationFile			text file containing annotation list {annotation; protein id list; protein symbol list}
	 * @return proteinToFrequencyMap	HashMap<String, Integer> {protein id: number of occurrence}
	 */
	public void computeFrequencyOfProteins(String annotationFile, String annotationCompanionFile, String outputFile){

		HashSet<String> motifsToTest = loadMotifsToTest(annotationCompanionFile);

		HashMap<String, Integer> proteinToFrequencyMap = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));

			String line = in.readLine(); // no header

			while(line!=null) {

				if(motifsToTest.contains(line.split("\\t")[0])) {

					String[] protein_ids = line.split("\\t")[2].split("\\|"); // idx[2] = protein (name) list 
					ArrayList<String> proteinList = checkProteinsInNetwork(protein_ids);

					if(proteinList.size() >= lowerBoundToSample && proteinList.size() <= upperBoundToSample) {
						/* For all proteins of a given annotation; if in list update number of occurrence, otherwise initialize */
						for(int i=0; i<proteinList.size(); i++) {
							if(proteinToFrequencyMap.containsKey(proteinList.get(i))) {
								proteinToFrequencyMap.replace(proteinList.get(i), proteinToFrequencyMap.get(proteinList.get(i)) + 1);
							} else {
								proteinToFrequencyMap.put(proteinList.get(i), 1);
							}
						}
					}
				}
				line = in.readLine();
			}			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		printProteinsToOccurence(outputFile, proteinToFrequencyMap);
	}

	private HashSet<String> loadMotifsToTest(String annotationCompanionFile){

		HashSet<String> motifSet = new HashSet<>();

		try {
			InputStream in = new FileInputStream(new File(annotationCompanionFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				motifSet.add(line.split("\t")[0]);
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}

	/**
	 * Output annotated proteins and their occurrences in a text file
	 * @param outputFile				text file to contain Protein : Occurrence
	 * @param proteinToOccurrenceMap	map of {protein: occurrence}
	 */
	private static void printProteinsToOccurence(String outputFile, HashMap<String, Integer> proteinToOccurrenceMap) {
		try {
			OutputStreamWriter out = new OutputStreamWriter(new FileOutputStream(new File(outputFile)));
			/* iterate through map */
			for(String protein: proteinToOccurrenceMap.keySet()) {
				out.write(protein + "\t" + proteinToOccurrenceMap.get(protein) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private ArrayList<String> checkProteinsInNetwork(String[] proteinList) {

		ArrayList<String> finalProteinList = new ArrayList<>();

		for(String prot: proteinList) {
			if(this.proteinSet.contains(prot)) {
				finalProteinList.add(prot);
			}
		}

		return finalProteinList;
	}

	public void combineProteinFrequencyData(String protFreqFilePrefix, int numFiles, String outputFile) {

		HashMap<String, Integer> proteinToOccurrenceMap = new HashMap<>();

		for(int i=0; i<numFiles; i++) {

			/* Load individual file protein frequencies and update proteinToOccurenceMap */
			String protFreqFile = protFreqFilePrefix + i;
			try {
				InputStream in = new FileInputStream(new File(protFreqFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine();

				while(line != null) {
					String[] col = line.split("\t");

					if(proteinToOccurrenceMap.containsKey(col[0])) {
						proteinToOccurrenceMap.put(col[0], proteinToOccurrenceMap.get(col[0]) + Integer.parseInt(col[1]));
					} else {
						proteinToOccurrenceMap.put(col[0], Integer.parseInt(col[1]));
					}
					line = input.readLine();
				}

				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

			printProteinsToOccurence(outputFile, proteinToOccurrenceMap);
		}
	}

}