package combineOutput;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class CombineLocalizationEnrichmentOutput {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/structures/localizations/";
		String enrichmentFilePrefix = "safeEnrichment_MotifBG_output_motif";
		String ouptputEnrichmentInformation = "safe/safeEnrichment_MotifBG_combinedData_motif";

		/* for each motif create a combined file - showing enrichment scores for significant domains (adjusted p-value < 0.05) */ 

		char[] dominantStructures = new char[] {'B', 'E', 'H', 'I', 'M', 'S', 'X'}; 
		for(int i=1; i<=31; i++) {

			HashMap<Character, StructureEnrichment> structures = new HashMap<>();

			for(char s: dominantStructures) {

				/* load enrichment information for given structure */
				
				String enrichmentFile = wd + "Motif" + i + "/" + enrichmentFilePrefix + i + "_" + s + ".tsv";
				File f = new File(enrichmentFile);
				if(f.exists() && !f.isDirectory()) { 
					structures.put(s, loadStructureEnrichment(enrichmentFile, s));
				}
			}
			/* print combined results */
			printEnrichmentResults(wd + ouptputEnrichmentInformation + i + ".tsv", structures, dominantStructures);
		}

	}


	public static StructureEnrichment loadStructureEnrichment(String inputFile, char s){

		StructureEnrichment struct;
		List<Double> significantScores = new ArrayList<>();
		List<Double> enrichmentScores = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line!=null) {

				String[] col = line.split("\t");
				significantScores.add(Double.parseDouble(col[7]));
				enrichmentScores.add(Double.parseDouble(col[5]));

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		struct = new StructureEnrichment(s, significantScores, enrichmentScores);
		return struct;
	}

	public static void printEnrichmentResults(String outputFile, HashMap<Character, StructureEnrichment> structures, char[] dominantStructures) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("SAFE\t");
			for(int i=1; i<=24; i++) {
				out.write(i + "\t");
			}
			out.write("\n");
			
			for(int i=0; i<dominantStructures.length; i++) {
				
				if(structures.containsKey(dominantStructures[i])) {
					
					out.write(dominantStructures[i] + "\t");
					
					double[] enrichmentScores = structures.get(dominantStructures[i]).getFormattedEnrichmentScores();
					for(int j=0; j<enrichmentScores.length; j++) {
						out.write(enrichmentScores[j] + "\t");
					}
					out.write("\n");
				}
			
				out.flush();

			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
