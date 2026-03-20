package rbpdb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class FormatPWMForTOMTOMquery {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/RBPDB/";
		
		// folder with position weighted matrice
		String pwmFolder = wd + "motifs/significantModules_pwm_";
		int[] motifs = new int[]{82,85, 94, 239};
		
		String outputFiile = wd + "LESMoNLocal_nwTPD2_queryMotifs_forTomtom.tsv";
		
		printFormattedQueryMotifsForTomTom(pwmFolder, motifs, outputFiile);
		
	}

	private static void printFormattedQueryMotifsForTomTom(String pwmPrefix, int[] motifs, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("MEME version 4\n\n" + "ALPHABET= ACGU\n\n" + 
					"Background letter frequencies\nA 0.25 C 0.25 G 0.25 U 0.25\n"); // header
			
			// for each motif file in folder
			for(Integer i : motifs) {
				System.out.print(i + ".");

				String pwmFile = pwmPrefix + i + ".tsv";
				/* Load the ppm (position probability matrix) for every family*/
				double[][] ppm = loadPPM(pwmFile);
				out.write("MOTIF" + i + "\n");
				out.write("letter-probability matrix: alength= 4 w= " + ppm[0].length + "\n");

				for(int k=0; k<ppm[0].length; k++) {
					for(int j=0; j<ppm.length; j++) {
						out.write(ppm[j][k] + "\t");
					}
					out.write("\n");
					out.flush();
				}
				out.write("\n");
			}
			out.write("\n"); // end signal

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Done");
	}
	
	private static double[][] loadPPM(String inputFile) {

		try {
			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			int lineCount = 0;

		    String[] col = line.split("\t");
			double[][] ppm = new double[4][col.length];
			
			while(line != null) {
				col = line.split("\t");

				for(int i=0; i<col.length; i++) {
					ppm[lineCount][i] = Double.parseDouble(col[i]);
				}

				line = input.readLine();
				lineCount++;
			}
			input.close();
			
			return ppm;
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
}
