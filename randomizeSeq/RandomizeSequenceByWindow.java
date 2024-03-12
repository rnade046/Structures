import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class RandomizeSequenceByWindow {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/Structures/hcm/";

		//		String fastaFile = wd + "input_files\\human_3UTRsequences.txt";
		//		String rnaIdListFile = wd + "motif_enumeration\\BiomaRt_MappingRefSeqIdsToGeneSymbol_corrNet.tsv";
		//		String randomFastaFile = wd + "input_files\\random_humanCellMap_3UTRsequences.txt";

		String fastaFile = wd + "corrNet2-400_3utr_w100cds.fasta";
		String randomFastaFile = wd + "corrNet2-400_3utr_w100cds_winShuffled_rand1.fasta";

		System.out.println("**Generating randomized fasta sequences**");
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(randomFastaFile)));
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));
			int countId = 0;

			System.out.println("Randomizing fasta");

			String line = in.readLine();
			while(line!=null) {

				if(countId%10 == 0) {
					System.out.print(countId + ".");
				} 
				if(countId%100 == 0) {
					System.out.println();
				}

				if(line.startsWith(">")) {
					out.write(line + "_Shuffled\n");
				} else {

					/* Get formatted sequence as String */
					String cds = randomizeSequence(line.substring(0,100));
					String utr = randomizeSequence(line.substring(100));

					out.write(cds + utr + "\n");
					countId++;
				}
				out.flush();
				line = in.readLine();
			}
			
			System.out.println("\nRandomized sequences: " + countId);
			in.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * For a given sequence, generate a it's random sequence using non overlapping window
	 * of 10 character 
	 * 
	 * @param seq	String - sequence to randomize
	 * @return randomSeq 	String - randomized sequence
	 */
	public static String randomizeSequence(String seq) {

		String randomSeq = "";

		/* iterate over sequence with a non overlaping window of size 10 */
		for(int i=0; i<seq.length(); i+=10) {

			String substring = "";

			/* substring of length 10, unless remainder sequences contains less then 10 character */ 
			if((seq.length() - i) >= 10) {
				substring = seq.substring(i, (i+10)); // returns substring of length 10
			} else {
				substring = seq.substring(i); // returns remainder of sequence (length < 10) 
			}

			if(substring.length() >= 2) {

				int max = substring.length();
				int min = 0;

				/* pairwise switch 1000 times two characters */ 
				for(int j=0; j<1000; j++) {
					String currentSubstring = "";

					/* initialize indexed of characters to switch, cannot be the same index or map to the same characters */
					int idx1, idx2 = 0;
					do {
						idx1 = (int) ((Math.random() * (max - min)) + min);
						idx2 = (int) ((Math.random() * (max - min)) + min);
					} while(idx1 == idx2);

					/* generate new substring, switching characters at selected indexes */
					for(int k=0; k<substring.length(); k++) {
						if(k == idx1) {
							currentSubstring += substring.charAt(idx2);
						} else if (k == idx2) {
							currentSubstring += substring.charAt(idx1);
						} else {
							currentSubstring += substring.charAt(k);
						}
					}
					substring = currentSubstring; // re-initialize substring
				}
			}

			randomSeq += substring; // add substring to seq  
		}
		return randomSeq;
	}
}
