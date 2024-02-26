package assessBayesPairing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class FilterAnnotationFile {

	public static void main(String[] args) {

		String annotationFile = args[0];
		String condition = args[2];

		double threshold = Double.parseDouble(args[1]);

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File("corrNetTop2-400_structureModules_" + condition + ".tsv")));

			String line = in.readLine(); // header
			out.write(line + "\n");

			line = in.readLine();
			while(line != null) {

				String[] col = line.split("\t"); // [0]=module, [1]=#Prots, [2]=listProteins
				out.write(col[0] + "\t");

				if(col.length > 2) {
					/* filter protein list based on score threshold */
					int protCount = 0;
					String proteinList="";
					for(String prot: col[2].split("\\|")) {

						String[] values = prot.split("\\_"); //[0]=proteinName, [1]=score
						if(Double.parseDouble(values[1]) >= threshold) {
							proteinList += prot + "|";
							protCount++;
						}
					}
					out.write(protCount + "\t" + proteinList + "\n");
					System.out.println("Filtered Mod " + col[0] + ":" + protCount + "/" + col[1]);
				} else {
					out.write("\n");
				}
				out.flush();
				line = in.readLine();
			}
			out.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
