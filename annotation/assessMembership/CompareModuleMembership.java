package assessMembership;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

public class CompareModuleMembership {

	public static void main(String[] args) {

		String fwdModulesMembershipFile="/Users/rnadeau2/Documents/Structures/Annotations/modules/membership/corrNet2-400_fwdSeq_bp4_proteinsByModule.tsv";
		String randModulesMembershipFile="/Users/rnadeau2/Documents/Structures/Annotations/modules/membership/corrNet2-400_seqRand1_bp4_proteinsByModule.tsv";
		
		String membershipOverlapFile = "/Users/rnadeau2/Documents/Structures/Annotations/modules/membership/membershipOverlap_seq_bp4.tsv";

		/* load module membership matrix */	
		HashMap<String, Module> modules1 = loadModules(fwdModulesMembershipFile);
		HashMap<String, Module> modules2 = loadModules(randModulesMembershipFile);

		/* compare membership */
		calculateMembershipOverlap(modules1, modules2, membershipOverlapFile);
	}

	public static HashMap<String, Module> loadModules(String membershipFile){

		HashMap<String, Module> moduleMap = new HashMap<>();

		/* load matrix */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(membershipFile))));

			String line = in.readLine(); // header
			String[] modules = Arrays.copyOfRange(line.split("\t"), 1, line.split("\t").length);

			List<List<String>> proteinSets = new ArrayList<>();
			for(int i=0; i<modules.length; i++) {
				proteinSets.add(new ArrayList<>());
			}
			line = in.readLine();

			while(line!=null) {

				String[] elements = line.split("\t");

				for(int i=1; i<elements.length; i++) {
					if(elements[i].equals("1")){
						proteinSets.get(i-1).add(elements[0]);
					}
				}
				line = in.readLine();
			}
			in.close();

			/* store as module objects */
			for(int i=0; i<modules.length; i++) {
				moduleMap.put(modules[i], new Module(modules[i], proteinSets.get(i)));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return moduleMap;
	}

	public static void calculateMembershipOverlap(HashMap<String, Module> modules1, HashMap<String, Module> modules2, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
				
			/* header */
			out.write("Module\tnProteins(fwd)\tnProteins(rand1)\tOverlap\n");
			
			for(Entry<String, Module> m : modules1.entrySet()) {
				if(modules2.containsKey(m.getKey())) {
					
					HashSet<String> proteins1 = new HashSet<>(m.getValue().getProteins());
					HashSet<String> proteins2 = new HashSet<>(modules2.get(m.getKey()).getProteins());
					
					int count = 0;
					for(String p : proteins1) {
						if(proteins2.contains(p)) {
							count++;
						}
					}
					
					double membership = count / (double) Math.min(proteins1.size(), proteins2.size());
					out.write(m.getKey() + "\t" + proteins1.size() + "\t" + proteins2.size() + "\t" + membership + "\n");
					out.flush();
				}
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
