package extra;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class CombineRandomizedAnnotations {

	public static void main(String[] args) {

		String inputFilePrefix = "/Users/rnadeau2/Documents/Structures/corrNetTop2-400_structureModules_rand";
		String outputFile = "/Users/rnadeau2/Documents/Structures/corrNetTop2-400_structureModules_randCombined_0.0.tsv";
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("Module\t#Prot\tProteinList\n");
			
			for(int i=1; i<=5; i++) {
				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFilePrefix + i + "_0.0.tsv"))));
				
				String line = in.readLine(); // header
				line = in.readLine();
				 
				while(line!=null) {
					
					String[] value = line.split("\t");
					
					out.write(value[0] + "_" + i + "\t");
					
					for(int j=1; j<value.length; j++) {
						out.write(value[j] + "\t");
					}
					out.write("\n");
					
					out.flush();
					line = in.readLine();
				}
				
				in.close();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		
		
		
		
		
	}

}
