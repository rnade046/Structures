package reliable;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;

import org.json.JSONObject;


public class AnalyseRELIABLEmodules {

	public static void main(String[] args) {

		String inputFile = "/Users/rnadeau2/Documents/Structures/BayesPairs/RELIABLE/RELIABLE.json";
		String outputFile = "/Users/rnadeau2/Documents/Structures/BayesPairs/RELIABLE/reliable-info.tsv";
		loadJson(inputFile, outputFile);
		
	}

	public static void loadJson(String inputFile, String outputFile) {

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			String line = in.readLine();
			
			out.write("module\tmodule_name\tPDBs\n");

			JSONObject rootObject = new JSONObject(line); // Parse the JSON to a JSONObject
			Iterator<String> keys = rootObject.keys();

			while(keys.hasNext()){

				/* load modules */
				String k = keys.next();
				
				JSONObject m = rootObject.getJSONObject(k); 
				
				JSONObject pdbs = m.getJSONObject("PDBs");
				String name = m.get("atlas_name").toString();
			
				out.write(k + "\t" + name + "\t");
				for(String pdb : pdbs.keySet()) {
					out.write(pdb + "|");
				}
				
				out.write("\n");
				out.flush();
			}
			in.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
