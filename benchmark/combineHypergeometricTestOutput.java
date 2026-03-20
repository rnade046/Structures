import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class combineHypergeometricTestOutput {

	public static void main(String[] args) {

		String inputPrefix = "/Users/rnadeau2/Documents/Structures/benchmark/i4/moduleInfo_output_clustSize3_module";
		String combinedDataFile = "/Users/rnadeau2/Documents/Structures/benchmark/output/corrNet2-400_benchmark_structures_score4_clustSize3_i4.tsv";

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(combinedDataFile)));
			
			/* header */
			out.write("Module\t");
			for(int i=1; i<=221; i++) {
				out.write("c" + i + "\t");
			}
			out.write("\n");
			
			/* for all modules 0 - 270 */
			for(int i=0; i <= 270; i++) {

				
				File f = new File(inputPrefix + i + ".tsv");
				
				List<Double> pvals = new ArrayList<Double>(130);
				if(f.exists()) {
					/* input significant p-values */
					pvals = loadPvalues(inputPrefix + i + ".tsv");
					
				}
				
				/* print values */
				if(!pvals.isEmpty()) {
					out.write("m" + i + "\t");
					for(int j=0; j<pvals.size(); j++) {
						out.write(pvals.get(j) + "\t");
					}
					out.write("\n");
					out.flush();	
				}
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public static List<Double> loadPvalues(String moduleInfoFile){

		List<Double> pvalues = new ArrayList<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(moduleInfoFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {
				double pval = Double.parseDouble(line.split("\t")[7]);
				pvalues.add(pval);

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return pvalues;
	}
}
