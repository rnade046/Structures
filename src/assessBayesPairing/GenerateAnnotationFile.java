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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import org.json.JSONArray;
import org.json.JSONObject;

public class GenerateAnnotationFile {

	public static void main(String[] args) {

		/* load protein list, but only modules with score that passes our set threshold */

		String wd = args[0];
		
		double threshold = Double.parseDouble(args[1]);
		String condition = args[2];
		
		String jsonIdxFile = "jsonIdxOfRefSeqIds_" + condition + ".tsv";
		String refSeqIdFile = "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String annotationFile = "corrNetTop2-400_structureModules_" + condition + "_" + threshold + ".tsv";

		/* determine mapping of protein - refSeqIds - .JSON files */ 
		List<Protein> proteinList = determineProteinMapping(refSeqIdFile, jsonIdxFile, condition);

		System.out.println("searching modules");
		for(Protein prot: proteinList) {
			System.out.println(prot);

			for(String jsonFileSuffix : prot.getFileIdMap().keySet()) {

				String refSeq = prot.getFileIdMap().get(jsonFileSuffix);
				String jsonFile = wd + jsonFileSuffix;
				System.out.println("file: " + jsonFile);
				File f = new File(jsonFile);
				if(f.exists()) {
					System.out.println(refSeq + " - checked");
					prot.updateMultipleModules(loadJson(jsonFile, refSeq, threshold));
				} else {
					System.out.println(refSeq + " - file not found");
					prot.addMissedId(refSeq); // add missed ID to list
				}
			}
			prot.summarizeModules();
		}
		/* create annotation file */
		generateAnnotationFile(proteinList, annotationFile);

	}

	public static List<Protein> determineProteinMapping(String refSeqToProteinFile, String refSeqToJSONFile, String condition) {

		List<Protein> proteinList = new ArrayList<>();

		System.out.println("loading json IdxMap");
		HashMap<String, String> jsonIdxMap = loadJSONIdx(refSeqToJSONFile);

		System.out.println("loading refSeqIds");
		HashMap<String, String[]> refSeqIdMap = loadRefSeqId(refSeqToProteinFile);

		System.out.println("Creating protein list:");
		for(Entry<String, String[]> entry: refSeqIdMap.entrySet()) {

			System.out.print(proteinList.size() + ".");

			if(proteinList.size()%50 == 0) {
				System.out.println();
			}

			Protein p = new Protein(entry.getKey());
			
			/* search refSeqIds */ 
			if(condition.equals("rand")) {
				for(String id: entry.getValue()) {
					p.updateJSONmapping(id + "_Shuffled", jsonIdxMap.get(id + "_Shuffled"));
				}
			} else { 
				for(String id: entry.getValue()) {
					p.updateJSONmapping(id, jsonIdxMap.get(id));
				}
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
				if(col.length > 1) {
					refSeqIdMap.put(col[0], col[1].split("\\|"));	
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

	public static List<Module> loadJson(String inputFile, String id, double threshold) {

		List<Module> modules = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();

			JSONObject rootObject = new JSONObject(line); // Parse the JSON to a JSONObject

			JSONObject a = rootObject.getJSONObject("all_hits").getJSONObject(id); // {270}

			Iterator<String> keys = a.keys();

			while(keys.hasNext()){

				/* load modules */
				String k = keys.next();
				JSONArray m = a.getJSONArray(k); 

				for(int i=0; i < m.length(); i++) {

					/* each instance of a module */
					JSONArray elements = m.getJSONArray(i);
					String seq = elements.get(0).toString();
					String[] pos = elements.get(1).toString().split("[,\\[\\]]");
					String[] pos2 = Arrays.copyOfRange(pos, 2, pos.length);
					double score = Double.parseDouble(elements.get(2).toString());

					if(score > threshold) { // positive score
						if(!pos2[0].isEmpty()) {
							if(Integer.parseInt(pos2[0]) > 100) { // ignore models in the CDS
								modules.add(new Module(Integer.parseInt(k), seq, score, pos2));
							}	
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

	public static void generateAnnotationFile(List<Protein> proteinList, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("Module\t#Prot\tProteinList\n");

			for(int i=0; i<=270; i++) {

				HashSet<String> protSet = new HashSet<>();

				for(Protein p : proteinList) {
					if(p.getModuleSummaryMap().containsKey(i)) {
						protSet.add(p.getProteinName());
					}
				}
				
				out.write(i + "\t" + protSet.size() + "\t");
				for(String prot: protSet) {
					out.write(prot + "|");
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
