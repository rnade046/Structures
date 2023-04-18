package sequences;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class FormatRNAfoldOutput {

	public static void main(String[] args) {

		String rnaFoldOutputFile = args[0];
		String formattedFile = args[1];

		/* format RNA fold output file to remove the score at the end of structure line */ 
		formatRNAfoldOutput(rnaFoldOutputFile, formattedFile);

	}

	public static void formatRNAfoldOutput(String inputFile, String outputFile) {

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			String line = in.readLine();
			int lineCount = 1;
			
			
			while(line!= null) {
				
				if(lineCount%3 == 0) {
					String structure = line.split("\\s+")[0];
					out.write(structure + "\n");
					
				} else {
					out.write(line + "\n");
				}
				
				out.flush();
				
				line = in.readLine();
				lineCount++;
			}
			
			in.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
