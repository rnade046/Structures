package localizations;

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

public class FormatHCMlocalizationForEnrichment {

	public static void main(String[] args) {
		
		String wd = "/Users/rnadeau2/Documents/Structures/";
		
		String proteinsInNetworkFile = wd + "enrichmentAnalysis/ioFiles/corrNetTop2-400_proteinsInNetwork_info.tsv";
		String preyFile = wd + "nwTPD/localization/preys-latest.txt"; 
		String annotationFile = wd + "enrichmentAnalysis/inputFiles/corrNetTop2-400_structureModules_fwd_4.0.tsv";
		String outputFile = wd + "nwTPD/localization/SAFE/safeEnrichmentInput_motif";
		String outputNMFFile = wd + "nwTPD/localization/NMF/nmfEnrichmentInput_motif";
		
		double threshold = 0.000902548470845385;
		String clusteringInfoFile = wd + "enrichmentAnalysis/fdr/bp4/corrNetTop2-400_nwTPD_fwd_bp4_structure_clusteringDetails.tsv";
		formatLocalizations(proteinsInNetworkFile, preyFile, clusteringInfoFile, threshold, annotationFile, outputFile, outputNMFFile);
		
	}


	public static void formatLocalizations(String proteinsInNetworkFile, String preyFile, String clusteringInfoFile, double threshold, String annotationFile, String outputFile, String outputNMFFile) {

		/* get proteins in network SAFE domain */
		System.out.println("** load prey info **");
		List<Protein> proteinList = loadPreyInfo(preyFile, loadProteinsInNetwork(proteinsInNetworkFile));
		
		/* determine output of each motif */
		
		/* load motifs */
		System.out.println("** load motifs to test **");
		HashSet<String> motifs = loadSignificantMotifs(clusteringInfoFile, threshold);
		
		/* for each motif get list of proteins */ 
		System.out.println("** get motif info **");
		proteinList = determineMotifsAssociatedToProteins(annotationFile, motifs, proteinList);
		
		/* print protein summary per motif */
		System.out.println("** print motifs **");
		printSafeDomains(outputFile, motifs, proteinList);
		printNMFregions(outputNMFFile, motifs, proteinList);
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

	/**
	 * Load representative motifs identified in the motif family process
	 *  
	 * @param motifFamilyFile	String - file containing all the representative motifs
	 * @return motifSet			Set<String> - set of representative motifs
	 */
	private static HashSet<String> loadSignificantMotifs(String motifFamilyFile, double threshold){

		HashSet<String> motifSet = new HashSet<>();
		try {

			InputStream in = new FileInputStream(new File(motifFamilyFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line!=null) {
				if(Double.parseDouble(line.split("\t")[3]) <= threshold) {
					motifSet.add(line.split("\t")[0]); // [0] = module number
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}	
	
	private static List<Protein> determineMotifsAssociatedToProteins(String annotationFile, HashSet<String> motifFamilies, List<Protein> proteinList){
		
		/* index protein list */
		HashMap<String, Integer> proteinIdxMap = new HashMap<>();
		for(int i=0; i<proteinList.size(); i++) {
			proteinIdxMap.put(proteinList.get(i).getName(), i);
		}
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			
			String line = in.readLine();
			while(line!=null) {
				
				String motif = line.split("\t")[0];
				if(motifFamilies.contains(motif)) {
					
;					String[] proteins = line.split("\t")[2].split("\\|");
					
					for(String p : proteins) {
						
						String pFormatted = p.split("\\_")[0];
						if(proteinIdxMap.containsKey(pFormatted)) {
							proteinList.get(proteinIdxMap.get(pFormatted)).addMotif(Integer.parseInt(motif));
						}
					}
					
				}
				
				line = in.readLine();
			}				
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return proteinList;
	}
	
	private static void printSafeDomains(String outputFile, HashSet<String> motifFamilies, List<Protein> proteinList) {
		
		for(String motif: motifFamilies) {
			
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + motif + ".tsv")));

				out.write("protein\tgroup\tSAFE\n");
				
				for(Protein prot : proteinList) {
					out.write(prot.getName() + "\t");
					
					if(prot.containsMotif(Integer.parseInt(motif))) {
						out.write("Motif\t");
					} else {
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
	
private static void printNMFregions(String outputFile, HashSet<String> motifFamilies, List<Protein> proteinList) {
		
		for(String motif: motifFamilies) {
			
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + motif + ".tsv")));

				out.write("protein\tgroup\tSAFE\n");
				
				for(Protein prot : proteinList) {
					out.write(prot.getName() + "\t");
					
					if(prot.containsMotif(Integer.parseInt(motif))) {
						out.write("Motif\t");
					} else {
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
}
