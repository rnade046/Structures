package goEnrichment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

public class FunctionalEnrichment {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		System.out.println("**Loading parameters file** \n");

		File directory = new File("ontologizer/"); 
		if (! directory.exists()){
			System.out.println("creating directory: ontologizer/");
			directory.mkdir();
		}

		/**
		 * ARGS[1] == 0 ; prepFunctionalEnrichment (b4 ONTOLOGIZER)
		 * ARGS[1] == 1 ; postFunctionalEnrichment (after ONTOLOGIZER; b4 REVIGO)
		 * ARGS[1] == 2 ; summarizeCellularComponents (after REVIGO)
		 */

		/* Load motif family info file to guide further analysis */
		String motifFamilyInfoFile = "corrNetTop2-400_clusteringDetails.tsv";

		System.out.println("**Loading representative motifs**");
		HashSet<String> motifSet = loadSignificantMotifs(motifFamilyInfoFile);
		System.out.println("Loaded motifs: " + motifSet.size());


		/* prepFunctionalEnrichment */
		if(args[0].equals("0")) {

			/* Generate file of all proteins in network for analysis - if it doesn't exist */ 
			String protAnnotationFreqFile = "ioFiles/corrNetTop2-400_proteinAnnotationFrequency.tsv";
			String backgroundFile = "ontologizer/corrNetTop2-400_background_annotatedProteins.txt";

			/* Create list of proteins in network from protein frequency file */
			System.out.println("**Format proteins in network**");
			formatProteinsInNetworkFile(protAnnotationFreqFile, backgroundFile);

			String annotatedProteinsPrefix = "ontologizer/annotatedProteinsByModule_";
			String annotationFile = args[2];

			String ontologizerIdxFile = "ontologizer/ontologizer-module-idx.tsv";
			/* Create 1 file per motif with the list of its annotated proteins for Ontologizer analysis */
			System.out.println("Format annotated proteins by motif:");
			formatAnnotatedProteinsByMotif(motifSet, annotationFile, annotatedProteinsPrefix, ontologizerIdxFile);
		}


		String tablePrefix = "ontologizer/table-annotatedProteinsByModule_";
		String outputPrefix = "ontologizer/significantGOterms_";

		/* postFunctionalEnrichment */
		if(args[0].equals("1")) {
			/* Load table files - get GO-terms associated to an adjusted p-value < 0.05 */

			System.out.println("**Extracting significant GO-terms from all proteins**");
			getSignificantGOterms(tablePrefix, outputPrefix, motifSet.size());

		}

		/* after REVIGO */
		if(args[0].equals("2")) {

			String revigoPrefix = "ontologizer/revigo_cc_";
			String cellularComponentSummaryFile = "ontologizer/cellularComponentSummary_matrix.tsv";

			/* Get significant cellular components */
			System.out.println("** Obtain significant cellular components - all proteins**");
			HashMap<String, String> ccInfoMap = getAllCellularComponents(revigoPrefix, motifSet.size());

			/* Summarize cellular components & print */
			System.out.println("**Summarize cellular compartment significance scores - all proteins**");

			List<List<GeneOntology>> significantCC = loadSignificantCellularComponents(tablePrefix, motifSet.size(), ccInfoMap);
			HashMap<String, Double[]> enrichmentMap = summarizeCellularComponentsEnrichmentScores(significantCC);
			printCellularComponentSummary(ccInfoMap, enrichmentMap, motifSet.size(), cellularComponentSummaryFile);

		}
	}

	private static void formatProteinsInNetworkFile(String proteinAnnotatedFreqFile, String proteinNetworkOutputFile) {

		try {

			InputStream in = new FileInputStream(new File(proteinAnnotatedFreqFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(proteinNetworkOutputFile)));

			String line = input.readLine();

			while(line!=null) {

				out.write(line.split("\t")[0] + "\n"); // protein name
				out.flush();
				line = input.readLine();
			}
			input.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Load representative motifs identified in the motif family process
	 *  
	 * @param motifFamilyFile	String - file containing all the representative motifs
	 * @return motifSet			Set<String> - set of representative motifs
	 */
	private static HashSet<String> loadSignificantMotifs(String motifFamilyFile){

		HashSet<String> motifSet = new HashSet<>();
		try {

			InputStream in = new FileInputStream(new File(motifFamilyFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line!=null) {
				motifSet.add(line.split("\t")[0]); // [0] = module number
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}	

	private static void formatAnnotatedProteinsByMotif(HashSet<String> repMotifSet, String extractedAnnotationFile, String annotatedProteinByMotifPrefix, String ontologizerIdxFile) {


		HashMap<String, Integer> moduleOrder = new HashMap<>();
		try {

			InputStream in = new FileInputStream(new File(extractedAnnotationFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // motif \t Prot1|Prot2|Prot3
			int motifCount = 1; 

			while(line!=null && motifCount<= repMotifSet.size()) {

				/* output necessary info when line in file corresponds to representative motif */
				if(repMotifSet.contains(line.split("\t")[0])) {

					System.out.println(motifCount);
					BufferedWriter out = new BufferedWriter(new FileWriter(new File(annotatedProteinByMotifPrefix + motifCount)));
					moduleOrder.put(line.split("\t")[0], motifCount);
					motifCount++;

					String[] protList = line.split("\t")[2].split("\\|");

					for(String prot: protList) {
						out.write(prot.split("\\_")[0] + "\n"); // protein name 
						out.flush();
					}
					out.close();
				}
				line = input.readLine(); // next line 
			}
			input.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
		printModuleOrder(moduleOrder, ontologizerIdxFile);

	}

	private static void printModuleOrder(HashMap<String,Integer> moduleOrder, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("Module\tOntologizer-Idx\n");
			for(Entry<String, Integer> module: moduleOrder.entrySet()) {
				out.write(module.getKey() + "\t" + module.getValue() + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void getSignificantGOterms(String tableInputFilePrefix, String outputPrefix, int numberOfFiles) {


		for(int i=1; i<= numberOfFiles; i++) {
			HashMap<String, String> goMap = new HashMap<>();
			try {

				InputStream in = new FileInputStream(new File(tableInputFilePrefix + i + "-Parent-Child-Union-Benjamini-Hochberg.txt"));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // header
				line = input.readLine(); 

				while(line!=null) {
					String[] col = line.split("\t"); // [0] = GO-term; [10] = adjusted p-val
					if(Double.parseDouble(col[10]) < 0.1) {
						goMap.put(col[0], col[10]);
					}

					line = input.readLine(); // next line 
				}
				input.close();

				if(!goMap.isEmpty()) {
					BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputPrefix + i)));

					for(Entry<String, String> e: goMap.entrySet()) {
						out.write(e.getKey() + "\t" + e.getValue() + "\n");
						out.flush();
					}

					out.close();
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	@SuppressWarnings("unused")
	private static void printSignificantGOterms(String outputPrefix, List<GeneOntology> significantGOList) {

		try {

			if(!significantGOList.isEmpty()) {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputPrefix)));

				for(GeneOntology go : significantGOList) {
					out.write(go.getName() + "\t" + go.getPval() + "\n");
					out.flush();
				}

				out.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static HashMap<String, String> getAllCellularComponents(String revigoPrefixFile, int numFamilies) {

		HashMap<String, String> ccInfo = new HashMap<>(); //GO:XX, Name

		for(int i=1; i<= numFamilies; i++) {
			System.out.print(i + ".");
			try {

				InputStream in = new FileInputStream(new File(revigoPrefixFile+ i));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // header
				for(int lineCount=1; lineCount<=5; lineCount ++) {
					line = input.readLine(); 
				}

				while(line!=null) {

					String[] col = line.split(","); // [0] = GO-term; [1] = name
					String term = col[0].split("\"")[1];
					String name = col[1].split("\"")[1];
					if(!ccInfo.containsKey(term)) {
						ccInfo.put(term, name);
					}
					line = input.readLine(); // next line 
				}
				input.close();

			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println("Done");
		return ccInfo;
	}		

	private static List<List<GeneOntology>> loadSignificantCellularComponents(String tableInputFilePrefix, int numFamilies, HashMap<String, String> ccInfo) {

		List<List<GeneOntology>> allSigngificantGOs = new ArrayList<>();

		for(int i=1; i<=numFamilies; i++) {
			List<GeneOntology> significantGOList = new ArrayList<>();

			try {

				InputStream in = new FileInputStream(new File(tableInputFilePrefix + i + "-Parent-Child-Union-Benjamini-Hochberg.txt"));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // header
				line = input.readLine(); 

				while(line!=null) {
					String[] col = line.split("\t"); // [0] = GO-term; [10] = adjusted p-val
					if(ccInfo.containsKey(col[0]) && Double.parseDouble(col[10]) < 0.05) {
						significantGOList.add(new GeneOntology(col));
					}

					line = input.readLine(); // next line 
				}
				input.close();

			} catch (IOException e) {
				e.printStackTrace();
			}

			allSigngificantGOs.add(significantGOList);
		}


		return allSigngificantGOs;
	}

	@SuppressWarnings("unused")
	private static HashMap<String, Double[]> summarizeCellularComponents(HashMap<String, String> ccInfoMap, String significantGOsPrefixFile, int numberFamily) {

		HashMap<String, Double[]> significantCCMap = new HashMap<>(); //GO:XXD, Double[family1, 2, .., x] (adjusted p-value)

		for(int i=1; i<=numberFamily; i++) {
			System.out.print(i + ".");
			try {

				InputStream in = new FileInputStream(new File(significantGOsPrefixFile+ i));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // header

				while(line!=null) {
					String[] col = line.split("\t"); // [0] = GO-term; [1] = adjusted p-val

					/* store info if significant GO:XX is a CC */
					if(ccInfoMap.containsKey(col[0])) {

						/* create entry if it's the first time this GO:XX is seen */ 
						if(!significantCCMap.containsKey(col[0])) {
							significantCCMap.put(col[0], new Double[numberFamily]);
						}
						/* update adjusted p-value at appropriate position */
						Double[] array = significantCCMap.get(col[0]);
						array[i-1] = Double.parseDouble(col[1]);

						significantCCMap.put(col[0], array);
					}
					line = input.readLine(); // next line 
				}
				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println("Done");
		return significantCCMap;
	}

	private static HashMap<String, Double[]> summarizeCellularComponentsEnrichmentScores(List<List<GeneOntology>> allSignificantCellularComponents) {

		HashMap<String, Double[]> significantCCMap = new HashMap<>(); //GO:XXD, Double[family1, 2, .., x] (adjusted p-value)

		for(int i=0; i<allSignificantCellularComponents.size(); i++) {
			System.out.print(i + ".");

			List<GeneOntology> currentMotif = allSignificantCellularComponents.get(i);

			for(GeneOntology go: currentMotif) {
				String goName = go.getName();

				/* create entry if it's the first time this GO:XX is seen */ 
				if(!significantCCMap.containsKey(goName)) {
					significantCCMap.put(goName, new Double[allSignificantCellularComponents.size()]);
				}
				/* update adjusted p-value at appropriate position */
				Double[] array = significantCCMap.get(goName);
				array[i] = go.getScore();

				significantCCMap.put(goName, array);
			}
		}
		System.out.println("Done");
		return significantCCMap;
	}

	private static void printCellularComponentSummary(HashMap<String, String> ccInfoMap, HashMap<String, Double[]> significantCCMap, int numFamilies, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			// header 
			out.write("GO-term\tGO-name\t#Groups\t");
			for(int i=1; i<= numFamilies; i++) {
				out.write(i + "\t");
			}
			out.write("\n");

			for(Entry<String, Double[]> goEntry: significantCCMap.entrySet()) {
				String go = goEntry.getKey();
				Double[] scores = goEntry.getValue();

				int count = 0;
				for(int i=0; i<scores.length; i++) {
					if(scores[i] != null) {
						count++;
					}
				}

				out.write(go + "\t" + ccInfoMap.get(go) + "\t" + count + "\t");
				for(int i=0; i<scores.length; i++) {
					out.write(scores[i] + "\t");
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
