package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class Scale {

	public static void scaleScores(String annotationFile, String scaledAnnotationFile) {

		/* determine current range (minimum and maximum) */

		Double[] rangeI = determineRange(annotationFile);
		Double[] rangeF = new Double[] {1.0, 10.0}; 

		/* format file */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(scaledAnnotationFile)));

			String line = in.readLine(); // header 
			out.write(line);
			
			line = in.readLine();

			while(line != null) {

				String[] values = line.split("\t");
				if(values.length > 2) {
					String[] proteins = values[2].split("\\|");

					out.write(values[0] + "\t" + values[1] + "\t");

					for(String prot: proteins) {
						String name = prot.split("\\_")[0];

						Double scoreI = Double.parseDouble(prot.split("\\_")[1]);
						
						/* scale score */
						Double scoreF = ((scoreI - rangeI[0])/(rangeI[1] - rangeI[0])) * (rangeF[1] - rangeF[0]) + rangeF[0];

						out.write(name + "_" + scoreF + "|");
					}
					out.write("\n");
					out.flush();
				}
				line = in.readLine();
			}
			out.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static Double[] determineRange(String annotationFile) {

		Double[] range = new Double[2];
		range[0] = Double.MAX_VALUE;
		range[1] = 0.0;

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));

			String line = in.readLine(); // header 
			line = in.readLine();

			while(line!=null) {

				/* Module	#Proteins	Protein1_score1|Protein2_score2|..|ProteinN_scoreN */
				if(line.split("\t").length > 2) {
					String[] values = line.split("\t")[2].split("\\|"); 

					for(String v: values) {
						double score = Double.parseDouble(v.split("\\_")[1]);

						/* check minimum */
						if(score < range[0]) {
							range[0] = score;
						}
						/* check maximum */
						if(score > range[1]) {
							range[1] = score;
						}
					}
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return range;
	}
}
