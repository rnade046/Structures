package checkSequenceStructure;

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
import java.util.regex.Matcher;

public class determineSequenceStructure {

	public static void main(String[] args) {


		String annotationFile = "corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";
		String refSeqIdFile = "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String motifsFile = "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";

		String indexFile = args[0];

		String structureFile = args[1];

		int motifToTest = Integer.parseInt(args[2]);

		/* index sequence files - ID = file name */
		HashMap<String, String> indexMap = getSequenceFilesByRefSeqId(indexFile);

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifsFile))));

			/* Read through significant motifs */

			String line = in.readLine(); // header 
			for(int i=0; i<motifToTest; i++) {
				line = in.readLine();
			}

			String motif = line.split("\t")[0];

			/* Set motifs as regex and restrict search for core proteins */
			Motif m = new Motif(motif, annotationFile, refSeqIdFile);

			System.out.println("* Motif info *");
			System.out.println("Loaded motif : " + m.getMotif());
			System.out.println("Regex pattern : " + m.getMotifRegexPattern() + "\n");
			
			/* load sequences - annotated bpRNA file - ignore first 100 bp - search for motif */ 
			for(String id: m.getRefSeqId().keySet()) {
				
				/* get the file name that contains RefSeqId corresponding sequence */
				if(indexMap.containsKey(id)) {
					String file = "annotatedStructures/" + indexMap.get(id) + ".st";

					/* search file for sequence info */
					searchStructureFilesForMotif(m, file, id);	
				}
				
				System.out.println("RefSeqId : " + id + "(number of occurrences: " + m.getOccurrences().size() + ")");
			}

			printMotifInfo(m, structureFile + motifToTest + ".tsv");
			line = in.readLine();

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static HashMap<String, String> getSequenceFilesByRefSeqId(String indexFile){

		HashMap<String, String> indexMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(indexFile))));

			String line = in.readLine();
			while(line!=null) {

				String[] entry = line.split("\t"); // [0] = File name, [1] = RefSeqId
				indexMap.put(entry[1],	entry[0]);

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return indexMap;
	}

	private static void searchStructureFilesForMotif(Motif m, String structureFile, String id) {

		List<Occurrence> motifOccurrences = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(structureFile))));


			String seq = "";
			String structure = "";

			for(int i=0; i<6; i++) {
				String line = in.readLine();

				if(i==3) {
					seq += line.subSequence(100, line.length());
				}

				if(i==5) {
					structure += line.substring(100, line.length());
				}
			}

			/* search sequence */ 
			Matcher matcher = m.getMotifRegexPattern().matcher(seq);	// match pattern to sequence 

			/* check for all instances of motif in sequence */
			while(matcher.find()) {

				/* obtain start position of motif */
				int pos = matcher.start(); 
				int end = pos+7;

				motifOccurrences.add(new Occurrence(m.getProtein(id), id, seq.substring(pos, end), structure.substring(pos, end),
						new int[] {pos, end}));
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		m.updateOccurrences(motifOccurrences);
	}

	private static void printMotifInfo(Motif m, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			//header 
			out.write("protein\tID\tsequenceMotifs\tStructureMotif\tDominantStructure\tStructureTypesCount\tPosition\n");

			for(Occurrence oc : m.getOccurrences()) {
				out.write(oc.getInfo() + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
