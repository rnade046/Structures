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

		String sequenceFile = "";
		String motifsFile = "";

		/* Read through significant motifs ; determine possible motif instances - map {Motif; Set<Instances>} */

		/* load sequence - annotated bpRNA file */ 


		/* Identify position of motif
		 * Identify what structures it's associated to
		 * Determine if it's split structure */



	}

	private HashMap<String, HashSet<String>> loadMotifsOfInterest(String motifsFile){
		HashMap<String, HashSet<String>> motifMap = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifsFile))));
			
				/* load motif family */
			
			
				/* determine set of possible motifs */

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifMap;
	}

}
