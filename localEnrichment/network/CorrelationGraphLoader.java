package network;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import graph.Interaction;

public class CorrelationGraphLoader {

	public static ArrayList<Interaction> loadGraphFromCorrelationNetwork(String inputRepository, String fastaFile, String mapProtToRefSeqIdsFile) {

		/* Load list of Proteins and their possible RefSeq IDs keeping only the IDs for which we have a sequence */
		HashMap<String, ArrayList<String>> mapProtToRefSeqIds = getRefSeqIdsInNetwork(mapProtToRefSeqIdsFile, fastaFile);

		/* Check number of interactions each protein is involved in; return Set of proteins to remove */
		 ArrayList<Interaction> confidentInteractions = formatInteractionList(inputRepository, mapProtToRefSeqIds);

		return confidentInteractions;
	}

	/**
	 * Get the list of refSeqIds associated to a given protein. RefSeqIds need to have a corresponding sequence in the fasta file
	 * 
	 * @param biomartMappingFile		String - file path output from BiomaRt outlining map of HGNC symbol to RefSeqId
	 * @param refSeqSet					String - Fasta file with sequences
	 * @return mapRefSeqIds				HashMap<String, ArrayList<String>> - map of protein as HGNC symbol = list of RefSeqIds
	 */
	private static HashMap<String, ArrayList<String>> getRefSeqIdsInNetwork(String biomartMappingFile, String fastaFile){

		HashMap<String, ArrayList<String>> mapRefSeqIds = new HashMap<>();

		HashSet<String> refSeqSet = generateRefSeqSet(fastaFile);

		try {
			InputStream in = new FileInputStream(new File(biomartMappingFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			while(line!=null) {

				String protein = line.split("\t")[0];
				String refSeqID = line.split("\t")[1];

				/* Add RefSeqId only if it has corresponding sequence in fasta file */
				if(refSeqSet.contains(refSeqID)) {

					/* Add refSeqId to existing list if protein has been seen b4, 
					 * or create new listing if protein hasn't been seen b4 */ 
					if(mapRefSeqIds.containsKey(protein)) {
						mapRefSeqIds.get(protein).add(refSeqID);
					} else {
						ArrayList<String> refSeqIdList = new ArrayList<>();
						refSeqIdList.add(refSeqID);
						mapRefSeqIds.put(protein, refSeqIdList);
					}
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return mapRefSeqIds;
	}

	/**
	 * Generate an index for the FASTA file; identifying the line number for the different refseq identifiers
	 * 
	 * @param fastaFile String - file path for the FASTA sequences
	 * @return indexOfFastaFile HashMap<String, Integer> - map of {refseqId : line count} 
	 */
	private static HashSet<String> generateRefSeqSet(String fastaFile){
		HashSet<String> refSeqSet = new HashSet<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(fastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				/* store the line index of the start of a new sequence */
				if(line.startsWith(">")) {
					String[] col = line.split("_|\\.");
					String refSeqId = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)

					refSeqSet.add(refSeqId);
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqSet;
	}

	private static ArrayList<Interaction> formatInteractionList(String networkFile, HashMap<String, ArrayList<String>> mapProtToRefSeqIds ){

		ArrayList<Interaction> interactionList = new ArrayList<>(); // list to contain formated interactions
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(networkFile))));

			String line = in.readLine(); // header
			line = in.readLine();
			
			while(line != null) {
				String[] col = line.split("\t");
				interactionList.add(new Interaction(col[0], col[1], mapProtToRefSeqIds.get(col[0]), mapProtToRefSeqIds.get(col[1]), (1-Double.parseDouble(col[2]))));
				
				line = in.readLine();
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return interactionList;
	}
}
