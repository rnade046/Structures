package network;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import graph.Interaction;

public class CorrelationGraphLoader {

	public static ArrayList<Interaction> loadGraphFromCorrelationNetwork(String inputRepository, String mapProtToRefSeqIdsFile) {

		/* Load list of Proteins and their possible RefSeq IDs keeping only the IDs for which we have a sequence */
		HashMap<String, ArrayList<String>> mapProtToRefSeqIds = getRefSeqIdsInNetwork(mapProtToRefSeqIdsFile);

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
	private static HashMap<String, ArrayList<String>> getRefSeqIdsInNetwork(String mappingFile){

		HashMap<String, ArrayList<String>> mapRefSeqIds = new HashMap<>();


		try {
			InputStream in = new FileInputStream(new File(mappingFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			while(line!=null) {

				String protein = line.split("\t")[0];
				
				if(line.split("\t").length > 1) {
					String[] refSeqID = line.split("\t")[1].split("\\|");
					mapRefSeqIds.put(protein, new ArrayList<String>(Arrays.asList(refSeqID)));
				} else {
					mapRefSeqIds.put(protein, new ArrayList<>());
				}

				line = input.readLine();
			}
			
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return mapRefSeqIds;
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
