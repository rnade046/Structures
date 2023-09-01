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
import java.util.List;

public class AssessProteinAnnotations {

	private int lowerBoundToSample;
	private int upperBoundToSample;

	private HashSet<String> proteinSet;

	public AssessProteinAnnotations(int _lowerBoundToSample, int _upperBoundToSample, HashSet<String> _proteinSet) {
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

		List<Protein> annotatedProteins = new ArrayList<>(); // info about annotated proteins 
		HashMap<String, Integer> proteinIndex = new HashMap<>(); // index of proteins in annotated protein list (initialized above)

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));

			String line = in.readLine(); // no header

			while(line!=null) {

				if(motifsToTest.contains(line.split("\\t")[0])) {

					String[] protein_ids = line.split("\\t")[2].split("\\|"); // idx[2] = protein (name) list + scores

					HashSet<String> proteinList = checkProteinsInNetwork(protein_ids);

					if(proteinList.size() >= lowerBoundToSample && proteinList.size() <= upperBoundToSample) {

						/* For all proteins of a given annotation; if in list update number of occurrence, otherwise initialize */
						for(int i=0; i<protein_ids.length; i++) {

							String prot = protein_ids[i].split("\\_")[0];

							if(proteinIndex.containsKey(prot)) {
								Protein currentProtein = annotatedProteins.get(i);

								currentProtein.updateFrequency(); // update protein count
								currentProtein.addScores(Double.parseDouble(protein_ids[i].split("\\_")[1])); // update scores
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

		printProteinsToOccurence(outputFile, annotatedProteins);
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
	private static void printProteinsToOccurence(String outputFile, List<Protein> proteinInfo) {
		try {
			OutputStreamWriter out = new OutputStreamWriter(new FileOutputStream(new File(outputFile)));
			/* iterate through map */
			for(Protein p : proteinInfo) {
				out.write(p.getProteinName() + "\t" + p.getAnnotationFrequency() + "\t");

				for(Double score: p.getScores()) {
					out.write(score + "|");
				}
				out.write("\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private HashSet<String> checkProteinsInNetwork(String[] proteinList) {

		HashSet<String> finalProteinList = new HashSet<>();

		for(String prot: proteinList) {

			String formattedProt = prot.split("\\_")[0];

			if(this.proteinSet.contains(formattedProt)) {
				finalProteinList.add(formattedProt);
			}
		}

		return finalProteinList;
	}

	@SuppressWarnings("unused")
	private void combineProteinFrequencyData(String protFreqFilePrefix, int numFiles, String outputFile) {

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

			//printProteinsToOccurence(outputFile, proteinToOccurrenceMap);
		}
	}

}