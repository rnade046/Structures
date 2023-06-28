package assessBayesPairing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class IndexJSONFiles {

	public static void main(String[] args) {
		
		String inputPrefix = args[0]; // "struct/corrNet2-400_3utr_w100cds_formattedStructure_"
		String outputFile = args[1];
		int start = Integer.parseInt(args[2]);
		int end = Integer.parseInt(args[3]);
		String condition = args[4];
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(int i=start; i<=end; i++) {
				
				if(i%50 == 0) {
					System.out.println();
				}
				System.out.print(i + "|");
				
				/* read .txt file */
				String inputFile = inputPrefix + i+".txt";
				
				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
				
				int refSeqCount = 0;
				String line = in.readLine();
				while(line !=null) {
					
					// >NM_0000000
					// >NM_0000000_Shuffled
					
					String refSeq = line.split(">")[1];
					String jsonPrefix = "bp_corrNet2-400_3utr_w100cds_" + condition + "_" + i + "_"+ refSeqCount+ ".json"; // format JSON file name 
					
					out.write(refSeq + "\t" + jsonPrefix + "\n");
					out.flush();
					
					// skip 3 lines
					for(int j=0; j<3; j++) {
						line = in.readLine();
					}
					refSeqCount++;
				}
				
			
				in.close();
			}
			
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}

}
