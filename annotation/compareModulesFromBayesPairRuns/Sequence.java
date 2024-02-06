package compareModulesFromBayesPairRuns;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import org.json.JSONArray;
import org.json.JSONObject;


public class Sequence {
	
	private String refSeqId;
	private String proteinName;
	
	private HashMap<String, Double> modulesMap1;
	private HashMap<String, Double> modulesMap2;
	
	private HashSet<String> newModules;
	private HashSet<String> lostModules;
	private HashSet<String> keptModules;
	
	public Sequence(String id, String prot, String jsonFile1, String jsonFile2) {
		
		this.refSeqId = id;
		this.proteinName = prot; 
		
		this.modulesMap1 = loadJson(jsonFile1);
		this.modulesMap2 = loadJson(jsonFile2);
		
		this.newModules = determineNewModules();
		this.lostModules = determineLostModules();
		this.keptModules = determineKeptModules();
	}

	private HashMap<String, Double> loadJson(String inputFile) {

		HashMap<String, Double> modulesMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();

			JSONObject rootObject = new JSONObject(line); // Parse the JSON to a JSONObject

			if(rootObject.getJSONObject("all_hits").has(this.refSeqId)){
				
				JSONObject a = rootObject.getJSONObject("all_hits").getJSONObject(this.refSeqId); // {270}

				Iterator<String> keys = a.keys();

				while(keys.hasNext()){

					/* load modules */
					String k = keys.next();
					JSONArray m = a.getJSONArray(k); 

					//System.out.println("k = " + k + " ; ");

					for(int i=0; i < m.length(); i++) {

						/* each instance of a module */
						JSONArray elements = m.getJSONArray(i);
						String[] pos = elements.get(1).toString().split("[,\\[\\]]");
						String[] pos2 = Arrays.copyOfRange(pos, 2, pos.length);
						double score = Double.parseDouble(elements.get(2).toString());

						if(score > 4) { // positive score
							if(!pos2[0].isEmpty()) {
								if(Integer.parseInt(pos2[0]) > 100) { // ignore models in the CDS
									modulesMap.put(k, score);
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
		return modulesMap;
	}
	
	private HashSet<String> determineNewModules(){
		
		HashSet<String> newModules = new HashSet<>();
		
		for(String module: this.modulesMap2.keySet()) {
			if(!this.modulesMap1.containsKey(module)) {
				newModules.add(module);
			}
		}
		
		return newModules;
	}
	
	private HashSet<String> determineLostModules(){
		
		HashSet<String> lostModules = new HashSet<>();
		
		for(String module: this.modulesMap1.keySet()) {
			if(!this.modulesMap2.containsKey(module)) {
				lostModules.add(module);
			}
		}
		return lostModules;
	}
		
	private HashSet<String> determineKeptModules(){
		
		HashSet<String> keptModules = new HashSet<>();
		
		for(String module: this.modulesMap1.keySet()) {
			if(this.modulesMap2.containsKey(module)) {
				keptModules.add(module);
			}
		}
		return keptModules;
	}
	
	public HashSet<String> getNewModules(){
		return this.newModules;
	}
	
	public HashSet<String> getLostModules(){
		return this.lostModules;
	}
	
	public HashSet<String> getKeptModules(){
		return this.keptModules;
	}
	
	public String getProteinName() {
		return this.proteinName;
	}
}
