package assessBayesPairing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.json.*;

public class assessBayesPairingOutput {

	public static void main(String[] args) {

		String wd = args[0];
		String condition = args[1];
		boolean shuffled = Boolean.parseBoolean(args[2]);

		String jsonIdxFile = "jsonIdxOfRefSeqIds_" + condition + ".tsv";
		String refSeqIdFile = "corrNetTop2-400_proteinsInNetwork_info.tsv";

		String outputFile = "ModuleScoreSummary-" + condition + ".tsv";
		String moduleRangeFile = "module-range-info-" + condition + ".tsv";
		String moduleInfoFile = "protein-module-structure-info-" + condition + ".tsv";
		String jsonInfoFile = "json-information-" + condition + ".tsv";

		String percentilesFile = "percentile-info-" + condition + ".tsv";

		/* determine mapping of protein - refSeqIds - .JSON files */ 
		List<Protein> proteinList = determineProteinMapping(refSeqIdFile, jsonIdxFile, shuffled);

		/* initialize map <JSON file : JSON entry> */
		HashMap<String, JSON> jsonMapping = initializeJSONmapping(jsonIdxFile);
		/* search for modules
		 * 
		 *  TO ADD : check if there's a protein without module-structure information ; maybe make a separate function to check this
		 *  */

		System.out.println("searching modules");
		int fileNotFound = 0; 
		for(Protein prot: proteinList) {
			System.out.println(prot);
			int count = 0;
			for(String jsonFileSuffix :prot.getFileIdMap().keySet()) {

				String jsonFile = wd + jsonFileSuffix;
				
				count++;
				System.out.println(prot.getFileIdMap().get(jsonFileSuffix) + " - checked: " + count + "/" + jsonMapping.size()); //
				System.out.println("jsonFileSuffix: " + jsonFileSuffix);
				
				if(jsonMapping.containsKey(jsonFileSuffix)) {
					jsonMapping.get(jsonFileSuffix).setProtein(prot.getProteinName());
				}
				File f = new File(jsonFile);
				if(f.exists()) {
					prot.updateMultipleModules(loadJson(jsonFile, prot.getFileIdMap().get(jsonFileSuffix)));

					if(jsonMapping.containsKey(jsonFileSuffix)) {
						jsonMapping.get(jsonFileSuffix).fileExists();
					}
				} else {
					System.out.println(prot.getFileIdMap().get(jsonFileSuffix) + " - file not found");
					prot.addMissedId(prot.getFileIdMap().get(jsonFileSuffix) ); // add missed ID to list
					fileNotFound++;
				}
			}
			prot.summarizeModules();
		}
		System.out.println("** Files not found : " + fileNotFound + " **");

		System.out.println("** determine score percentiles **");
		CalculatePercentileRanks.determinePercentiles(proteinList, percentilesFile);

		List<List<Double>> moduleSummary = combineModuleCounts(proteinList);

		System.out.println("** print module summary **");
		System.out.println("modules : " + moduleSummary.size());
		printCombinedResults(outputFile, moduleSummary);

		System.out.println("** protein info summary **");
		assessModuleStructureInfo(proteinList, moduleInfoFile);

		System.out.println("** print module range info **");
		printModuleRange(moduleRangeFile, moduleSummary, proteinList);

		/* print JSON file information */
		System.out.println("** print json info file **");
		printJSONinfo(jsonInfoFile, jsonMapping);
	}

	public static List<Protein> determineProteinMapping(String refSeqToProteinFile, String refSeqToJSONFile, boolean shuffled) {

		List<Protein> proteinList = new ArrayList<>();

		System.out.println("loading json IdxMap");
		HashMap<String, String> jsonIdxMap = loadJSONIdx(refSeqToJSONFile);

		System.out.println("loading refSeqIds");
		HashMap<String, String[]> refSeqIdMap = loadRefSeqId(refSeqToProteinFile, shuffled);

		System.out.println("Creating protein list:");
		for(Entry<String, String[]> entry: refSeqIdMap.entrySet()) {

			System.out.print(proteinList.size() + ".");

			if(proteinList.size()%50 == 0) {
				System.out.println();
			}

			Protein p = new Protein(entry.getKey());
			for(String id: entry.getValue()) {
				p.updateJSONmapping(id, jsonIdxMap.get(id));
			}

			proteinList.add(p);
		}
		return proteinList;
	}

	/** load map of RefSeqIds to protein names */
	public static HashMap<String, String[]> loadRefSeqId(String refSeqToProteinFile, boolean shuffled){

		HashMap<String, String[]> refSeqIdMap = new HashMap<>(); // Protein = {Id1, Id2, ..., IdN}
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(refSeqToProteinFile)));

			String line = in.readLine();
			while(line!=null) {
				String[] col = line.split("\t");
				if(col.length > 1) {

					String[] ids = col[1].split("\\|");
					if(shuffled) {
						String[] idShuffled = new String[ids.length];
						for(int i=0; i<ids.length; i++) {
							idShuffled[i] = ids[i] + "_Shuffled";
						}

						refSeqIdMap.put(col[0], idShuffled);
					} else {
						refSeqIdMap.put(col[0], ids);	
					}


				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqIdMap;
	}

	public static HashMap<String, String> loadJSONIdx(String refSeqToJSONFile){

		HashMap<String, String> jsonIdxMap = new HashMap<>(); // RefSeqId = JSON File 
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(refSeqToJSONFile)));

			String line = in.readLine();
			while(line!=null) {
				String[] col = line.split("\t");

				jsonIdxMap.put(col[0], col[1]);

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return jsonIdxMap;
	}

	public static List<List<Double>> combineModuleCounts(List<Protein> proteinList) {

		List<List<Double>> modules = new ArrayList<>(); 
		for(int i=0; i<=270; i++) { // 0-270 modules
			modules.add(new ArrayList<Double>());
		}

		for(Protein prot : proteinList) {

			for(Entry<Integer, Double> entry: prot.getModuleSummaryMap().entrySet()) {
				//System.out.println("module = " + entry.getKey());
				//System.out.println("value = " + entry.getValue());
				modules.get(entry.getKey()).add(entry.getValue()); // key = module, value = score
			}
		}
		return modules;
	}

	/* load map of JSON files to refSeqIds */
	public static void printCombinedResults(String outputFile, List<List<Double>> modules) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("Module\tS>0\tS>0.5\tS>1\tS>1.5\tS>2\tS>2.5\tS>3\tS>3.5\tS>4\n");

			for(int i=0; i<modules.size(); i++) {
				List<Double> moduleScores = modules.get(i);

				int[] count = new int[9]; 	// 0-4, interval of 0.5
				// [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]

				for(double score: moduleScores) {

					int index = 0;

					if(score > 4) {
						index = 8; 
					} else if(score > 3.5) {
						index = 7;
					} else if(score > 3) {
						index = 6;
					} else if(score > 2.5) {
						index = 5;
					} else if(score > 2) {
						index = 4;
					} else if(score > 1.5) {
						index = 3;
					} else if(score > 1) {
						index = 2;
					} else if(score > 0.5) {
						index = 1;
					}

					for(int j=0; j<=index; j++) {
						count[j] += 1;
					}
				}

				out.write(i +"\t");

				for(int j=0; j<count.length; j++) {
					out.write(count[j] +"\t");
				}

				out.write("\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
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

						if(score > 0) { // positive score
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

	public static void assessModuleCount(List<Module> modules, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(Module m : modules) {
				out.write(m.getID() + "\t" + m.getSequence() + "\t" + m.getScore() + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static void printModuleRange(String outputFile, List<List<Double>> modules, List<Protein> protList) {

		double modMin = Double.MAX_VALUE;
		double modMax = 0;

		double protMin = Double.MAX_VALUE;
		double protMax = 0;
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("Module\tMinScore\tMaxScore\n");
			for(int i=0; i<modules.size(); i++) {

				List<Double> mScores = modules.get(i);
				if(mScores.size() > 0) {

					double currentMin = Collections.min(mScores);
					double currentMax = Collections.max(mScores);

					if(currentMin < modMin) {
						modMin = currentMin;
					}

					if(currentMax > modMax) {
						modMax = currentMax;
					}

					out.write(i + "\t" + currentMin + "\t" + currentMax + "\n");
				}
				out.flush();
			}

			out.write("\n\n\n\n");
			out.write("Protein\tMinScore\tMaxScore\n");

			for(Protein p : protList) {
				/* get scores */
				ArrayList<Double> scores = new ArrayList<>(p.getModuleSummaryMap().values());

				if(!scores.isEmpty()) {
					double currentMin = Collections.min(scores);
					double currentMax = Collections.max(scores);

					if(currentMin < protMin) {
						protMin = currentMin;
					}

					if(currentMax > protMax) {
						protMax = currentMax;
					}

					out.write(p.getProteinName() + "\t" + currentMin + "\t" + currentMax + "\n");
					out.flush();
				}
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		System.out.println("Module | Min : " + modMin + " | Max : " + modMax);
		System.out.println("Protein | Min : " + protMin + " | Max : " + protMax);

	}


	public static void assessModuleStructureInfo(List<Protein> protList, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("Protein\t#RefSeqIDs\t#MissedIds\tListMissedIds\n");

			for(Protein p: protList) {
				out.write(p.getProteinName() + "\t" + p.getFileIdMap().size() + "\t" + p.getMissedIDs().size() + "\t");

				for(String id: p.getMissedIDs()) {
					out.write(id + "|");
				}
				out.write("\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static HashMap<String, JSON> initializeJSONmapping(String jsonIndexFile){

		HashMap<String, JSON> jsonMapping = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(jsonIndexFile)));

			String line = in.readLine();
			while(line!=null) {

				String[] col = line.split("\t");
				jsonMapping.put(col[1], new JSON(col[1], col[0]));

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return jsonMapping;
	}

	public static void printJSONinfo(String outputFile, HashMap<String, JSON> jsonMapping) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("JSONFile\tRefSeqId\tProtein\tSeqSize\tExists\n");

			for(JSON json : jsonMapping.values()) {

				String[] values = json.getEntry();
				for(int i=0; i<values.length; i++) {
					out.write(values[i] + "\t");
				}
				out.write("\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
