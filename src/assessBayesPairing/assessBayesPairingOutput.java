package assessBayesPairing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.json.*;

public class assessBayesPairingOutput {

	public static void main(String[] args) {

		String[] files = new String[] {"cds_structure", "cds", "default", "sample1000", "structure_size1000", "structure"} ; 

		for(String f : files) {
			String jsonFile = "/Users/rnadeau2/Documents/Structures/BayesPairs/output/" + f + ".json";
			String outputFile = "/Users/rnadeau2/Documents/Structures/BayesPairs/"+ f+".tsv";

			List<Module> modules = loadJson(jsonFile);
			assessModuleCount(modules, outputFile);

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

					if(score > 0) {
						modules.add(new Module(Integer.parseInt(k), seq, score, pos));
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
