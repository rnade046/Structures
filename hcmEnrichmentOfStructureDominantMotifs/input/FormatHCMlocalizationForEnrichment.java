package input;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

public class FormatHCMlocalizationForEnrichment {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";

		String proteinsInNetworkFile = wd + "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String preyFile = wd + "localization/preys-latest.txt"; 
		String motifsInterestFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv"; 

		String motifStructurePrefix = wd + "structures/motifStructures/motifStructures_";

		String outputFile = "/safeEnrichmentInput_motif";
		String outputNMFFile = "/nmfEnrichmentInput_motif";
		
		String outputFile2 = "/safeEnrichmentInput_MotifBG_motif";
		String outputNMFFile2 = "/nmfEnrichmentInput_MotifBG_motif";
		
		


		/* get proteins in network SAFE domain */
		System.out.println("** load prey info **");
		List<Protein> proteinList = loadPreyInfo(preyFile, loadProteinsInNetwork(proteinsInNetworkFile));

		/* load motifs {String motif = int motifIndex} */
		System.out.println("** load motifs to test **");
		HashMap<String, Integer> motifs = getMotifsToTest(motifsInterestFile);

		/* for each motif get list of proteins */ 
		System.out.println("** get motif info **");

		for(int m: motifs.values()) {

			File f = new File(wd + "structures/localizations/Motif" + m + "/"); 

			if (!f.exists()) { 
				f.mkdirs();
			} 

			/* load motif structure information */  
			MotifStructures structures = getMotifStructures(m, motifStructurePrefix);

			/* print protein summary per motif */
			System.out.println("** print motifs **");
			printSafeDomains(f + outputFile, m, proteinList, structures.getStructures());
			printNMFregions(f + outputNMFFile, m, proteinList, structures.getStructures());
			
			printSafeDomainsMotifOnly(f+ outputFile2, m, proteinList, structures.getStructures());
			printNMFregionsMotifOnly(f + outputNMFFile2, m, proteinList, structures.getStructures());
		}


	}

	private static HashSet<String> loadProteinsInNetwork(String proteinsInNetworkFile){

		HashSet<String> proteinSet = new HashSet<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinsInNetworkFile))));

			String line = in.readLine();
			while(line !=null) {

				proteinSet.add(line.split("\t")[0]);
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinSet;
	}

	private static List<Protein> loadPreyInfo(String preyFile, HashSet<String> proteinSet){ 

		List<Protein> proteinList = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(preyFile))));

			String line = in.readLine(); // header 
			line = in.readLine();

			while(line!=null) {

				String[] col = line.split("\t");

				if(proteinSet.contains(col[0])) {

					int safe = 25;
					int nmf = 19;

					if(!col[3].equals("-")) {
						safe = Integer.parseInt(col[3]);
					}

					if(!col[1].equals("19")) {
						nmf = Integer.parseInt(col[1]);
					}

					proteinList.add(new Protein(col[0], safe, nmf));
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinList;
	}

	private static HashMap<String, Integer> getMotifsToTest(String motifsInterestFile){

		HashMap<String, Integer> motifSet = new HashMap<>();

		try {
			InputStream in = new FileInputStream(new File(motifsInterestFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			int count = 1;

			while(line != null) {

				motifSet.put(line.split("\t")[0], count);
				line = input.readLine();
				count++;
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifSet;
	}

	private static MotifStructures getMotifStructures(int motifId, String motifStructureFilePrefix) {

		MotifStructures structures = new MotifStructures();

		try {
			InputStream in = new FileInputStream(new File(motifStructureFilePrefix + motifId + ".tsv"));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line != null) {

				String[] col = line.split("\t");
				structures.updateStructures(col[4].charAt(0), col[0]);

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return structures;
	}

	private static void printSafeDomains(String outputFile, int motif, List<Protein> proteinList, HashMap<Character, HashSet<String>> structures) {

		for(Entry<Character, HashSet<String>> s: structures.entrySet()) {

			char dStructure = s.getKey();
			HashSet<String> proteinsWithStructure = s.getValue();

			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + motif + "_" + dStructure + ".tsv")));

				out.write("protein\tgroup\tSAFE\n");

				for(Protein prot : proteinList) {


					out.write(prot.getName() + "\t");

					if(proteinsWithStructure.contains(prot.getName())) {
						out.write("Motif\t");
					}else {
						out.write("NotMotif\t");
					}

					out.write(prot.getSafeDomain() + "\n");
					out.flush();
				}

				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

	}

	private static void printSafeDomainsMotifOnly(String outputFile, int motif, List<Protein> proteinList, HashMap<Character, HashSet<String>> structures) {


		HashSet<String> proteinSet = new HashSet<>();
		for(HashSet<String> proteins: structures.values()) {
			proteinSet.addAll(proteins);
		}

		for(Entry<Character, HashSet<String>> s: structures.entrySet()) {
			char dStructure = s.getKey();
			HashSet<String> proteinsWithStructure = s.getValue();

			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + motif + "_" + dStructure + ".tsv")));

				out.write("protein\tgroup\tSAFE\n");

				for(Protein prot : proteinList) {

					if(proteinSet.contains(prot.getName())) {

						out.write(prot.getName() + "\t");

						if(proteinsWithStructure.contains(prot.getName())) {
							out.write("Motif\t");
						}else {
							out.write("NotMotif\t");
						}

						out.write(prot.getSafeDomain() + "\n");
						out.flush();
					}
				}

				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

		}

	}


	private static void printNMFregions(String outputFile, int motif, List<Protein> proteinList, HashMap<Character, HashSet<String>> structures) {

		for(Entry<Character, HashSet<String>> s: structures.entrySet()) {

			char dStructure = s.getKey();
			HashSet<String> proteinsWithStructure = s.getValue();

			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + motif + "_" + dStructure +".tsv")));

				out.write("protein\tgroup\tSAFE\n");

				for(Protein prot : proteinList) {
					out.write(prot.getName() + "\t");

					if(proteinsWithStructure.contains(prot.getName())) {
						out.write("Motif\t");
					}else {
						out.write("NotMotif\t");
					}

					out.write(prot.getNMFregion() + "\n");
					out.flush();
				}

				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

	private static void printNMFregionsMotifOnly(String outputFile, int motif, List<Protein> proteinList, HashMap<Character, HashSet<String>> structures) {

		HashSet<String> proteinSet = new HashSet<>();
		for(HashSet<String> proteins: structures.values()) {
			proteinSet.addAll(proteins);
		}

		for(Entry<Character, HashSet<String>> s: structures.entrySet()) {

			char dStructure = s.getKey();
			HashSet<String> proteinsWithStructure = s.getValue();

			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + motif + "_" + dStructure +".tsv")));

				out.write("protein\tgroup\tSAFE\n");

				for(Protein prot : proteinList) {

					if(proteinSet.contains(prot.getName())) {

						out.write(prot.getName() + "\t");

						if(proteinsWithStructure.contains(prot.getName())) {
							out.write("Motif\t");
						}else {
							out.write("NotMotif\t");
						}
						out.write(prot.getNMFregion() + "\n");
						out.flush();
					}
				}
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}
}
