package sequenceLength;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class assessSignificantModules {

	public static void main(String[] args) {

		String modulesDetailsFile = "/Users/rnadeau2/Documents/Structures/enrichmentAnalysis/fdr/bp4/corrNetTop2-400_nwTPD_fwd_bp4_structure_clusteringDetails.tsv";
		String annotationFile = "/Users/rnadeau2/Documents/Structures/enrichmentAnalysis/inputFiles/corrNetTop2-400_structureModules_fwd_4.0.tsv";

		String proteinInfoFile = "/Users/rnadeau2/Documents/Structures/enrichmentAnalysis/ioFiles/corrNetTop2-400_proteinsInNetwork_info.tsv";
		String fastaFile = "/Users/rnadeau2/Documents/Structures/hcm/corrNet2-400_3utr_w100cds.fasta";

		String idxFile = "/Users/rnadeau2/Documents/Structures/nwTPD/corrNet2-400_nwTPD_fwd_bp4_idxOfFasta.tsv";
		String fastaFilePrefix = "/Users/rnadeau2/Documents/Structures/nwTPD/fasta/nwTPD_fwd_bp4_sequence_";

		String proteinInfoFileFinal = "/Users/rnadeau2/Documents/Structures/nwTPD/corrNet2-400_nwTPD_fwd_bp4_proteinInfo.tsv";

		/* get list of significant modules */
		System.out.println("** load significant modules **");
		HashSet<Integer> modules = getSignificantModules(modulesDetailsFile, 0.000902548470845385);

		/* get list {module = list of proteins} */
		System.out.println("** load proteins of interest **");
		HashSet<String> proteinSet = getProteinsAnnotatedByModules(modules, annotationFile);

		/* get list {protein = refSeqIds} */
		System.out.println("** load refSeqIds **");
		List<Protein> proteins = getRefSeqIds(proteinInfoFile, proteinSet);

		/* obtain sequence lengths of refSeqIds -- update to print relevant sequences info for BayesPairing */
		System.out.println("** assess sequence lengths **");
		getFastaInfo(proteins, fastaFile, idxFile, fastaFilePrefix);

		System.out.println("** print information **");
		printProteinInformation(proteins, proteinInfoFileFinal);

	}

	public static HashSet<Integer> getSignificantModules(String inputFile, double threshold){

		HashSet<Integer> modules = new HashSet<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine(); // header

			while(line!=null) {
				String[] col = line.split("\t");

				if(Double.parseDouble(col[3]) <= threshold) {
					modules.add(Integer.parseInt(col[0]));
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return modules;
	}

	public static HashSet<String> getProteinsAnnotatedByModules(HashSet<Integer> modules, String annotationFile){

		HashSet<String> proteins = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));

			String line = in.readLine(); // header
			line = in.readLine();
			while(line!=null) {
				String[] col = line.split("\t");

				if(modules.contains(Integer.parseInt(col[0]))) {
					String[] protElements = col[2].split("\\|");
					for(String p : protElements) {
						proteins.add(p.split("\\_")[0]);
					}
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteins;
	}

	public static List<Protein> getRefSeqIds(String proteinInfoFile, HashSet<String> proteinSet){

		List<Protein> proteins = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfoFile))));

			String line = in.readLine(); // header

			while(line!=null) {
				String[] col = line.split("\t");

				if(proteinSet.contains(col[0])) {
					if(col.length>1) {
						proteins.add(new Protein(col[0], col[1].split("\\|")));
					}
				}

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteins;
	}

	public static void getFastaInfo(List<Protein> proteins, String fastaFile, String indexMappingFile, String fastaPrefixFile) {
		int protCount = 1; 
		int idFile = 0;

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(indexMappingFile)));
			out.write("idx\tID\tLength\tProtein\n");
			
			for(Protein p : proteins) {

				System.out.print(protCount+ ".");
				if(protCount%50==0) {
					System.out.println();
				}

				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));
				
				String currentID = "";
				String line = in.readLine(); 
				
				boolean loadSeq = false;

				while(line != null) { 

					if(loadSeq) {

						// update protein object with sequence length
						p.setSeqLength(currentID, line.length());

						// print index info
						out.write(idFile + "\t" + currentID + "\t" + line.length() + "\t" + p.getProteinName() + "\n");
						out.flush();

						// print FASTA
						printFastaFileForBayesPairing(fastaPrefixFile + idFile + ".fasta", currentID, line);

						loadSeq = false;
					}

					if(line.startsWith(">")) {

						//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
						currentID = line.split("\\>")[1];

						if(p.containsId(currentID)) {
							loadSeq = true;
							idFile++;
						}
					}
					line = in.readLine();
				}

				in.close();

				protCount++;
			}
			System.out.println();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static void printFastaFileForBayesPairing(String fastaFile, String id, String sequence) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(fastaFile)));

			out.write(">" + id + "\n"); // header 
			out.write(sequence + "\n"); // sequence

			out.flush();
			out.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void printProteinInformation(List<Protein> proteins, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("Protein\t#IDs\tList-IDs\tminLength\tmaxLength\n");

			for(Protein p : proteins) {

				out.write(p.getProteinName() +"\t"+ p.getNumberOfIDs() +"\t" + p.getIdInfo() + "\t" + p.getLengthLimits()[0] + "\t" + p.getLengthLimits()[1] + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
