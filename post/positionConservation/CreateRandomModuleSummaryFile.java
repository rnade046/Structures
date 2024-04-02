package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

public class CreateRandomModuleSummaryFile {

	public static void main(String[] args) {

		String wd = args[0];
		String fwdSignificantModuleFile = wd + args[1];
		String randJSONidxFile = wd + args[2];

		String randSignificantModuleFile = wd + "rand_correspondingModuleSummaryFile.tsv";

		/* map refSeqId to JSON indexes */
		HashMap<String, String> jsonIdxMap = loadJSONIdx(randJSONidxFile);

		/* create corresponding module info file from FWD significant module summary */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fwdSignificantModuleFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(randSignificantModuleFile)));
			
			/* header */
			String line = in.readLine();
			out.write(line);
			
			line = in.readLine();
			while(line!=null) {
				
				String[] info = line.split("\t");
				out.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[3] + "\t" + info[4] + "\t");
				
				for(String protein : info[5].split("\\|")) {
					
					out.write(protein.split("\\_",2)[0] + "_["); // protein name
					
					String[] ids = protein.split("\\_",2)[1].substring(1, protein.split("\\_", 2)[1].length()-1).split(",");
					for(String i: ids) { // refSeqId=jsonPath

						String refSeqId = i.split("=")[0] + "_Shuffled";
						if(jsonIdxMap.containsKey(refSeqId)) {
							out.write(refSeqId + "=" + jsonIdxMap.get(refSeqId)+ ",");
						}
					}
					out.write("]|");
				}
				out.write("\n");
				out.flush();
				line = in.readLine();
			}
			out.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
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

}
