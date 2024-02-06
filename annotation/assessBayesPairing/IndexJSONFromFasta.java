package assessBayesPairing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class IndexJSONFromFasta {

	public static void main(String[] args) {

		String fastaDirectory = args[0];
		String jsonPrefixFile = args[1];
		String jsonIdxFile = args[2];

		/* obtain all FASTA files */
		File[] listOfFiles = new File(fastaDirectory).listFiles();

		/* output {refSeqId : jsonFile name} */
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(jsonIdxFile)));

			/* for each FASTA - obtain */ 
			for(File f: listOfFiles) {

				/* obtain # id (remove initial zeroes) */
				String[] name = f.getName().split("[\\_\\.]");
				int id = Integer.parseInt(name[name.length-2]);
				System.out.println(id);
				/* obtain refSeq Id */
				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(f)));

				String line = in.readLine();

				// >NM_0000000
				// >NM_0000000_Shuffled
				String refSeq = line.split(">")[1];
				String jsonPrefix = jsonPrefixFile + id + ".json"; // format JSON file name 

				out.write(refSeq + "\t" + jsonPrefix + "\t");
				out.flush();

				in.close();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
