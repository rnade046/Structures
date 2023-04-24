package checkSequenceStructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

public class determineSequenceStructure {

	public static void main(String[] args) {


		String annotationFile = "";
		String refSeqIdFile = "";
		String motifsFile = "";
		
		String indexFile = "";

		String structureFile = "";

		/* index sequence files - ID = file name */
		HashMap<String, String> indexMap = getSequenceFilesByRefSeqId(indexFile);
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifsFile))));

			/* Read through significant motifs */

			String line = in.readLine(); // header 
			line = in.readLine();

			while(line !=null) {

				String motif = line.split("\t")[0];

				/* Set motifs as regex and restrict search for core proteins */
				Motif m = new Motif(motif, annotationFile, refSeqIdFile);

				/* load sequences - annotated bpRNA file - ignore first 100 bp - search for motif */ 
				for(String id: m.getRefSeqId().keySet()) {
					
					/* get the file name that contains RefSeqId corresponding sequence */
					String file = indexMap.get(id);
					
					/* search file for sequence info */
					
				}

				line = in.readLine();
			}

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
	
	private void searchStructureFilesForMotif(Motif m, String structureFile) {
		
		HashMap<String, String> indexMap = new HashMap<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(structureFile))));
			
			String line = in.readLine();
			while(line!=null) {
				
				/* load sequence */
				
				/* load structure */
			
				line = in.readLine();
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
