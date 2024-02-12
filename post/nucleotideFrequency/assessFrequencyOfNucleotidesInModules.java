package nucleotideFrequency;

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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;

import org.json.JSONArray;
import org.json.JSONObject;


public class assessFrequencyOfNucleotidesInModules {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		/* params for file paths?? ie. JSON & results need to be update*/
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));	

		String wd = params.getProperty("working_directory");
		String jsonPath = args[1];
		/* input files */
		String significantModulesFile = wd + "significantModules_summary.tsv";

		/* output */
		Files.createDirectories(Paths.get(wd + "nuclFreq/"));

		String moduleSequencesFilePrefix = wd + "nuclFreq/significantModules_seqInstances_"; // instances of sequences
		String modulePPMFilePrefix = wd + "nuclFreq/significantModules_pwm_"; // PWM

		/* get significant modules */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(significantModulesFile))));		/* search each module for sequence and structure frequency */

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				/* obtain module info */
				String[] info = line.split("\t");
				String module = info[0];
				System.out.println("testing module: " + module + " | #prots : " + info[4]);
				/* initialize list of sequences for protein */ 
				ArrayList<String> sequenceList = new ArrayList<String>();

				for(String protein : info[5].split("\\|")) {

					System.out.println(protein.split("\\",2)[0] + " : ");

					String[] ids = protein.split("\\_",2)[1].substring(1, protein.split("\\_", 2)[1].length()-1).split(",");

					Double score = 0.0;
					String sequence = "";

					for(String i: ids) { // refSeqId=jsonPath

						String[] identifiers = i.split("=");

						if(identifiers.length > 1) {
							File f = new File(jsonPath + identifiers[1]);
							if(f.exists()) {
								Object[] values = searchJsonForMaxModuleScore(jsonPath + identifiers[1], identifiers[0], module, score, sequence);

								if(values[0] != null) {
									score = Double.valueOf(values[0].toString());
									sequence = values[1].toString();
								}
							} 
						}
					}
					/* update list with module sequence that had >score for the given protein */
					if(!sequence.isBlank()) {
						sequenceList.add(sequence); 
						System.out.print("(" + score + ")");
					}
				}

				double[][] ppm = calculatePPM(sequenceList, sequenceList.get(0).length());

				printPPM(ppm, modulePPMFilePrefix + module + ".tsv");
				printMotifs(sequenceList, moduleSequencesFilePrefix + module + ".tsv");

				System.out.println();
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static Object[] searchJsonForMaxModuleScore(String inputFile, String id, String module, Double maxScore, String sequence) {

		Object[] r = new Object[2];

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();

			JSONObject rootObject = new JSONObject(line); // Parse the JSON to a JSONObject

			if(rootObject.getJSONObject("all_hits").has(id)){

				JSONObject a = rootObject.getJSONObject("all_hits").getJSONObject(id); // {270}
				System.out.print(id + ": ");

				if(a.has(module)) {

					JSONArray m = a.getJSONArray(module); 
					for(int i=0; i < m.length(); i++) {

						/* each instance of a module */
						JSONArray elements = m.getJSONArray(i);
						String seq = elements.get(0).toString();
						String[] pos = elements.get(1).toString().split("[,\\[\\]]");
						String[] pos2 = Arrays.copyOfRange(pos, 2, pos.length);
						double score = Double.parseDouble(elements.get(2).toString());

						if(score > 4) { // positive score
							if(!pos2[0].isEmpty()) {
								if(Integer.parseInt(pos2[0]) > 100) { // ignore models in the CDS

									if(score > maxScore) {
										r[0] = score;
										r[1] = seq;

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
		return r;
	}

	private static double[][] calculatePPM(ArrayList<String> motifInstances, int motifLength){

		double[][] pfm = new double[4][motifLength];
		double[][] ppm = new double[4][motifLength];

		/* compute position frequency matrix */
		for(String motif: motifInstances) {
			for(int i=0; i<motifLength; i++) {
				switch(motif.charAt(i)) {
				case 'A': 
					pfm[0][i] += 1;
					break;
				case 'C': 
					pfm[1][i] += 1;
					break;
				case 'G':
					pfm[2][i] += 1;
					break;
				case 'T':
					pfm[3][i] += 1;
					break;
				}
			}
		}

		/* convert position frequency matrix to position probability matrix */
		for(int i=0; i<4; i++) {
			for(int j=0; j<motifLength; j++) {
				ppm[i][j] = pfm[i][j]/ (double) motifInstances.size();
			}
		}
		return ppm;
	}

	private static void printPPM(double[][] ppm, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<ppm.length; i++) {
				for(int j=0;j<ppm[i].length; j++) {
					out.write(ppm[i][j] + "\t");
				}
				out.write("\n");
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printMotifs(ArrayList<String> motifInstances, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(String motif: motifInstances) {
				out.write(motif + "\n");
				out.flush();
			}


			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
