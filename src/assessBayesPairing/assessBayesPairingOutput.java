package assessBayesPairing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.json.*;

public class assessBayesPairingOutput {

	public static void main(String[] args) {

		String jsonIdxFile = "jsonIdxOfRefSeqIds.tsv";
		String refSeqIdFile = "corrNetTop2-400_proteinsInNetwork_info.tsv";

		String outputFile = "ModuleScoreSummary.tsv";

		/* determine mapping of protein - refSeqIds - .JSON files */ 
		List<Protein> proteinList = determineProteinMapping(refSeqIdFile, jsonIdxFile);

		/* search for modules */
		for(Protein prot: proteinList) {

			for(String jsonFileSuffix :prot.getFileIdMap().keySet()) {
				String jsonFile = "" + jsonFileSuffix;
				prot.updateMultipleModules(loadJson(jsonFile));
			}
			prot.summarizeModules();
		}
		List<List<Double>> moduleSummary = combineModuleCounts(proteinList);

		printCombinedResults(outputFile, moduleSummary);
	}

	public static List<Protein> determineProteinMapping(String refSeqToProteinFile, String refSeqToJSONFile) {

		List<Protein> proteinList = new ArrayList<>();

		HashMap<String, String> jsonIdxMap = loadJSONIdx(refSeqToJSONFile);
		HashMap<String, String[]> refSeqIdMap = loadRefSeqId(refSeqToProteinFile);

		for(Entry<String, String[]> entry: refSeqIdMap.entrySet()) {

			Protein p = new Protein(entry.getKey());
			for(String id: entry.getValue()) {
				p.updateJSONmapping(id, jsonIdxMap.get(id));
			}

			proteinList.add(p);
		}
		return proteinList;
	}

	/** load map of RefSeqIds to protein names */
	public static HashMap<String, String[]> loadRefSeqId(String refSeqToProteinFile){

		HashMap<String, String[]> refSeqIdMap = new HashMap<>(); // Protein = {Id1, Id2, ..., IdN}
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(refSeqToProteinFile)));

			String line = in.readLine();
			while(line!=null) {
				String[] col = line.split("\t");
				refSeqIdMap.put(col[0], col[1].split("\\|"));

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
		for(int i=0; i<196; i++) { // 0-196 modules
			modules.add(new ArrayList<Double>());
		}

		for(Protein prot : proteinList) {
			for(Entry<Integer, Double> entry: prot.getModuleSummaryMap().entrySet()) {
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
					out.write(j +"\t");
				}
				
				out.write("\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static List<Module> loadJson(String inputFile) {

		List<Module> modules = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();

			JSONObject rootObject = new JSONObject(line); // Parse the JSON to a JSONObject

			JSONObject a = rootObject.getJSONObject("all_hits").getJSONObject("input_seq"); // {270}

			Iterator<String> keys = a.keys();

			while(keys.hasNext()){

				/* load modules */
				String k = keys.next();
				JSONArray m = a.getJSONArray(k); 

				System.out.println("\nk = " + k + " ; ");

				for(int i=0; i < m.length(); i++) {

					/* each instance of a module */
					JSONArray elements = m.getJSONArray(i);
					String seq = elements.get(0).toString();
					String[] pos = elements.get(1).toString().split("[,\\[\\]]");
					double score = Double.parseDouble(elements.get(2).toString());

					if(score > 0) { // positive score
						if(Integer.parseInt(pos[0]) > 100) { // ignore models in the CDS
							modules.add(new Module(Integer.parseInt(k), seq, score, pos));
						}
					}
					System.out.print(i + " ");
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
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

}
