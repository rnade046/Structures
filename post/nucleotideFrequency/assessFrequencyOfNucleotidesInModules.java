package nucleotideFrequency;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.json.JSONArray;
import org.json.JSONObject;

import assessBayesPairing.Module;

public class assessFrequencyOfNucleotidesInModules {

	public static void main(String[] args) {
		
		/* input files */
		String annotationFile = "";
		String proteinInfoFile = ""; // map Protein - RefSeqIDs
		String jsonIdxFile = ""; // map RefSeqIds - BayesPair file
		
		String significantModulesFile = "";
		
		/* output */
		String sequenceFreqPrefixFile = ""; // PWM
		String structureFreqPrefixFile = ""; // instances of structures
		
		/* get significant modules */
		HashSet<String> significantModules = loadSignificantMotifs(significantModulesFile, 0.000902548470845385);
		
		/* set up modules {protein list, refSeqId list, .JSON file paths} */
		HashMap<String, Module> modules = initializeModules(significantModules, annotationFile, proteinInfoFile, jsonIdxFile);
		
		/* search each module for sequence and structure frequency */
	}
	
	private static HashSet<String> loadSignificantMotifs(String motifFamilyFile, double threshold){

		HashSet<String> motifSet = new HashSet<>();
		try {

			InputStream in = new FileInputStream(new File(motifFamilyFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line!=null) {
				if(Double.parseDouble(line.split("\t")[3]) <= threshold) {
					motifSet.add(line.split("\t")[0]); // [0] = module number
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}	
	
	public static HashMap<String, Module> initializeModules(HashSet<String> significantModules, String annotationFile, String proteinInfoFile, String jsonIdxFile){
		
		HashMap<String, String[]> idMap = loadProteinToIdMap(proteinInfoFile);
		HashMap<String, String> jsonIdxMap = loadJSONidxFile(jsonIdxFile);
		
		HashMap<String, Module> modules = new HashMap<String, Module>();
		
		/* load annotation file - keeping only info relevant to significant modules */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			String line = in.readLine();

			while(line!=null) {
				
				String[] col = line.split("\t");
				
				if(significantModules.contains(col[0])) {
					modules.put(col[0], new Module(col[2].split("\\|"), idMap, jsonIdxMap));
				}
				line = in.readLine();
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return modules;
	}
	
	public static HashMap<String, String[]> loadProteinToIdMap(String infoFile){
		
		HashMap<String, String[]> idMap = new HashMap<String, String[]>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(infoFile))));
			String line = in.readLine();

			while(line!=null) {
				
				String[] col = line.split("\t");
				
				if(col.length > 1) {
					idMap.put(col[0], col[1].split("\\|")); // [0] = protein, [1] = list of IDs
				}
				line = in.readLine();
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return idMap;
	}
	
	public static HashMap<String, String> loadJSONidxFile(String jsonIdxFile){
		
		HashMap<String, String> jsonIdxMap = new HashMap<String, String>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(jsonIdxFile))));
			String line = in.readLine();

			while(line!=null) {
				
				String[] col = line.split("\t");
				
				jsonIdxMap.put(col[0], col[1]); // [0] = id, [1] = JSON path
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return jsonIdxMap;
	}
	
	public static void searchForModuleInstances(HashMap<String, Module> modules, String sequenceFreqPrefixFile) {
		 
		for(Entry<String, Module> entry : modules.entrySet()) {
			
			ArrayList<String> sequenceInstances = new ArrayList<>();
			ArrayList<String> structureInstances = new ArrayList<>();			
			
			/* load JSON files individually - get entries corresponding to module */
		
			/* add entry to list */
			
			
			double[][] ppm = calculatePPM(sequenceInstances, sequenceInstances.get(0).length());
			printPPM(ppm, sequenceFreqPrefixFile + entry.getKey() + ".tsv");
			
			printMotifs(structureInstances, structureInstances + entry.getKey() + ".tsv");
		}
		

	}
	
	public static List<Module> loadJson(String inputFile, String id) {

		List<Module> modules = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();

			JSONObject rootObject = new JSONObject(line); // Parse the JSON to a JSONObject

			if(rootObject.getJSONObject("all_hits").has(id)){
				
				JSONObject a = rootObject.getJSONObject("all_hits").getJSONObject(id); // {270}

				Iterator<String> keys = a.keys();

				while(keys.hasNext()){

					/* load modules */
					String k = keys.next();
					JSONArray m = a.getJSONArray(k); 

					//System.out.println("k = " + k + " ; ");

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
									modules.add(new Module(Integer.parseInt(k), seq, score, pos2));
								}	
							}

						}
						System.out.print(i + " ");
					}
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		System.out.println("modules: " + modules.size());
		return modules;
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
