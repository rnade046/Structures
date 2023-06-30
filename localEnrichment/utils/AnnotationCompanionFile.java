package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

public class AnnotationCompanionFile {

	/**
	 * Load annotation file; check proteins associated to annotations 1-by-1; 
	 * if the actual number of proteins lies between the lower and upper bound, 
	 * print the motif, the number of originally listed proteins and the actual number of proteins (in this network)
	 * 
	 * @param annotationFile
	 * @param companionFile
	 * @param proteinSet
	 * @param numFiles
	 * @param lowerBound
	 * @param upperBound
	 */
	public static void determineAnnotatedProteinsInNetwork(String annotationFile, String companionFile, HashSet<String> proteinSet, int lowerBound, int upperBound) {

		try {
			InputStream in = new FileInputStream(new File(annotationFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(companionFile)));

			String line = input.readLine(); // no header

			while(line!= null) {

				String[] col = line.split("\t");
				String motif = col[0];
				int protCount = 0;
				String[] proteins = new String[0];
				
				if(col.length > 2) {
					proteins = col[2].split("\\|");

					/* Count proteins in network */

					for(String protein: proteins) {
						if(proteinSet.contains(protein)) {
							protCount++;
						}
					}
				}
				/* print motif info if it respects upper and lower bound */
				if(protCount >= lowerBound && protCount <= upperBound) {
					out.write(motif + "\t" + proteins.length + "\t" + protCount + "\n");
					out.flush();
				}

				line = input.readLine();
			}
			input.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
