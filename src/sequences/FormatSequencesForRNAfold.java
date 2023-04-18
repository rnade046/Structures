package sequences;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;

public class FormatSequencesForRNAfold {


	public static void main(String[] args) { 

		// files
		String utrSequencesFile = "/Users/rnadeau2/Documents/Structures/hcm/hg38_3utr_sequence.fasta.txt";
		String cdsSequencesFile = "/Users/rnadeau2/Documents/Structures/hcm/hg38_cds_sequence.fasta.txt";

		String genesInNetworkFile = "/Users/rnadeau2/Documents/Structures/hcm/MappingRefSeqIdsToGeneSymbol_corrNet2-400.tsv";
		
		String formatedSequencesFile = "/Users/rnadeau2/Documents/Structures/hcm/fasta/corrNet2-400_3utr_w100cds_";

		// load refSeqIds associated to proteins in network
		HashSet<String> idSet = loadIds(genesInNetworkFile);

		/* iterate over sequence IDs - get UTR and CDS sequence - format output for RNAfold */ 
		int count = 0;
		int fileCount = 0;
		for(String id: idSet) {
			
			System.out.print(count + ".");
			
			if(count%50 == 0) {
				System.out.println("\n");
			}
			
			String utr = getSequenceForGivenId(id, utrSequencesFile);
			String cds = getSequenceForGivenId(id, cdsSequencesFile);
			
			if(!utr.isEmpty() && !cds.isEmpty()) {
				printFormattedSequence(formatedSequencesFile + fileCount + ".fasta", id, utr, cds);
			}
			count++;
			
			if(count%14 == 0) {
				fileCount++; 
			}
		}

	}

	public static HashSet<String> loadIds(String file){

		HashSet<String> idSet = new HashSet<>();
		HashSet<String> geneSet = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(file))));

			String line = in.readLine(); 

			while(line != null) { 

				idSet.add(line.split("\t")[1]);
				geneSet.add(line.split("\t")[0]);

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Genes with sequence: " + geneSet.size());

		return idSet;
	}

	public static String getSequenceForGivenId(String id, String fastaFile) {

		String seq = "";

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));

			String line = in.readLine(); 
			boolean loadSeq = false;

			while(line != null) { 
				
				if(loadSeq & !line.startsWith(">")) {
					seq += line;
				}

				if(line.startsWith(">")) {
					
					loadSeq = false;

					/* if sequence has been loaded - exit loop and return sequence  */
					if(!seq.isEmpty()) {
						break;
					}

					//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
					String currentID = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];

					if(currentID.equals(id)) {
						loadSeq = true;
					}
					
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return seq;
	}
	
	public static void printFormattedSequence(String outputFile, String id, String utr, String cds) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile), true));
			
			out.write(">" + id + "\n" + cds.substring(cds.length()-100, cds.length()) + utr + "\n");
			
			out.flush();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
