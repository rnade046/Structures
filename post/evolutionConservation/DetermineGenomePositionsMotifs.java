package evolutionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;

import org.json.JSONArray;
import org.json.JSONObject;

public class DetermineGenomePositionsMotifs {

	public static void main (String[] args) {

		String wd = args[0];
		String jsonPath = args[1];
		
		String fastaFile = wd + args[2];
		String significantModulesFile = wd + args[3];

		String bedFile = wd + "genomeBED/genomePositions_module";
		
		/* create folder for bed file*/
		File directory = new File(wd + "genomeBED/");
		if(!directory.exists()) {
			directory.mkdir();
		}
		
		/* load 3UTR chromosomes positions Map<RefSeqId = Chr:strt-end> */
		HashMap<String, String> utrPositionMap = mapGenomicPositionsToRefSeqIds(fastaFile);

		/* iterate over significant modules */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(significantModulesFile))));		

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {

				/* obtain module info */
				String[] info = line.split("\t");
				String module = info[0];

				BufferedWriter out = new BufferedWriter(new FileWriter(new File(bedFile + module + ".bed")));
				out.write("Track - Module : " + module + "\n");
				
				System.out.println("testing module: " + module + " | #prots : " + info[4]);

				/* Determine genomic position of significant motif with highest score associated to each protein */
				for(String protein : info[5].split("\\|")) {

					System.out.println(protein.split("\\_",2)[0] + " : ");
					String[] ids = protein.split("\\_",2)[1].substring(1, protein.split("\\_", 2)[1].length()-1).split(",");

					/* initialize module for protein to compare against */
					Module utrMod = new Module(0, 0.0);
					String utrFinalIdentifiers = "";

					for(String i: ids) { // refSeqId=jsonPath

						String[] identifiers = i.split("=");

						/* search for module with highest score within each sequence corresponding to refSeqId */
						if(identifiers.length > 1) {
							File f = new File(jsonPath + identifiers[1]);
							if(f.exists()) {

								/* obtain information for current refSeqId*/
								Module utrCurrentMod = searchJsonForMaxModuleScoreInUTR(jsonPath + identifiers[1], identifiers[0], module, utrMod.getScore());

								/* update information if values array is not empty (i.e. corresponds to modules with > score) */
								if(utrCurrentMod != null) {
									utrMod = utrCurrentMod;
									utrFinalIdentifiers = identifiers[0]; // refSeqId
								}
							} 
						}
					} 

					/* append position of current module to bed file */
					if(utrMod.getScore() > 0) {
						System.out.println("\nutr > " + utrFinalIdentifiers + " | score: " + utrMod.getScore());
						
						String loci = utrPositionMap.get(utrFinalIdentifiers);
						out.write(determineMotifGenomicPosition(loci, utrMod) + "\n");
						out.flush();
						
					}
				}
				System.out.println();
				line = in.readLine();
				
				out.close();
			}
			in.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static HashMap<String, String> mapGenomicPositionsToRefSeqIds(String fastaFile){
		
		HashMap<String, String> positionMap = new HashMap<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));
			
			String line = in.readLine();
			while(line!=null) {
				
				/* check if sequence is associated to a protein in network */
				if(line.startsWith(">")) {
					
					//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
					String id = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];
					String pos = getUTRgenomicPosition(line);
					
					positionMap.put(id, pos);
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return positionMap;
	}

	private static String getUTRgenomicPosition(String header) {

		String chromosome = header.split(" ")[1].split("=")[1].split(":")[0];
		chromosome = chromosome+":";
		int startPos =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[0]);
		int endPos =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[1]);

		String genomicPosition = chromosome + startPos +"-"+endPos;

		return genomicPosition;
	}

	public static Module searchJsonForMaxModuleScoreInUTR(String inputFile, String id, String module, Double maxScore) {

		Module mod = null;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();

			JSONObject rootObject = new JSONObject(line); // Parse the JSON to a JSONObject
			int seqLength = rootObject.getJSONArray("input").get(0).toString().length();

			if(rootObject.getJSONObject("all_hits").has(id)){

				JSONObject a = rootObject.getJSONObject("all_hits").getJSONObject(id); // {270}
				//				System.out.print(id + ": ");

				if(a.has(module)) {

					JSONArray m = a.getJSONArray(module); 
					for(int i=0; i < m.length(); i++) {

						/* each instance of a module */
						JSONArray elements = m.getJSONArray(i);
						//String seq = elements.get(0).toString();
						String[] pos = elements.get(1).toString().split("[,\\[\\]]");
						String[] pos2 = Arrays.copyOfRange(pos, 2, pos.length);
						double score = Double.parseDouble(elements.get(2).toString());

						if(score > 4) { // positive score
							if(!pos2[0].isEmpty()) {
								if(Integer.parseInt(pos2[0]) > 100) { // ignore models in the CDS

									if(score > maxScore) {
										mod = new Module(seqLength-100, score, pos2, true);
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
		return mod;
	}
	
	private static String determineMotifGenomicPosition(String utrPosition, Module mod) {
		
		/* select module position at random */
		int[] positions = mod.getPositions();
		int randomIdx = (int) Math.random() * positions.length;
		
		/* obtain start position of utr */
		String[] loci = utrPosition.split(":");
		int start = Integer.parseInt(loci[1].split("-")[0]) + positions[randomIdx]; 
		
		String modulePosition = loci[0] +"\t"+ start + "\t" + (start+1);
		return modulePosition;
	}
}
