package reliable;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.Set;

import org.json.JSONObject;


public class AnalyseRELIABLEmodules {

	public static void main(String[] args) {

		String inputFile = "/Users/rnadeau2/Documents/Structures/BayesPairs/RELIABLE/RELIABLE.json";
		String outputFile = "/Users/rnadeau2/Documents/Structures/BayesPairs/RELIABLE/reliable-info2.tsv";
		loadJson(inputFile, outputFile);
	}

	public static void loadJson(String inputFile, String outputFile) {

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("Module\t3dAtlasID\tSiblings\tmGraph-nodes\tmGraph-edges\tPDBs\n");

			JSONObject rootObject = new JSONObject(in.readLine()); // Parse the JSON to a JSONObject
			Iterator<String> keys = rootObject.keys();

			while(keys.hasNext()){

				/* load modules */
				String k = keys.next();

				JSONObject m = rootObject.getJSONObject(k); 
				
				Set<String> pdbs = m.getJSONObject("PDBs").keySet();
				String name = m.get("atlas_name").toString();
				String sibs = m.getJSONArray("siblings").toString();
				
				JSONObject masterGraph = m.getJSONObject("master_graph");
				String nodes = masterGraph.get("nodes").toString();
				String edges = masterGraph.get("edges").toString();
				
				out.write(k + "\t" + name + "\t" + sibs.substring(1, sibs.length()-1) + "\t" + nodes.substring(1, nodes.length()-1) +
						"\t" + edges.substring(1, edges.length()-1) +"\t" + String.join(",", pdbs));
				
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
