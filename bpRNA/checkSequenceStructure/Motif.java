package checkSequenceStructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

public class Motif {

	private String motif; 
	private Pattern regexMotif;
	private HashMap<String, String> refSeqIdMap; // {ID = protein} 

	public Motif(String m, String annotationFile, String refSeqIdFile) {

		this.motif = m;
		this.regexMotif = Pattern.compile(formatMotifWithRegularExpression(m));
		this.refSeqIdMap = determineRefSeqIdsAssociatedToCoreProteinsAnnotatedByMotif(annotationFile, refSeqIdFile);
	}

	/* get functions */
	public String getMotif() {
		return this.motif;
	}
	
	public Pattern getMotifRegexPattern() {
		return this.regexMotif;
	}
	
	public HashMap<String, String> getRefSeqId(){
		return this.refSeqIdMap;
	}

	/**
	 * Take input motif and translate to corresponding regular expression. 
	 * 
	 * @param motif				String - initial motif
	 * @param characterMap		HashMap<Character, String> - conversion map with regular expressions
	 * @return formattedMotif	String - motif formatted with regular expression
	 */
	private String formatMotifWithRegularExpression(String motif) {

		String formattedMotif = "";

		/* set character map */
		HashMap<Character, String> characterMap = setCharacterMapForRegularExpression();

		/* determine regular expression */
		for(int i = 0; i < motif.length(); i++){
			formattedMotif += characterMap.get(motif.charAt(i));
		}

		return formattedMotif;
	}

	/**
	 * Set character map - hard coded with alphabet used throughout analysis 
	 * 
	 * Values follow known IUPAC regular expression e.g. R = "[AG]" 
	 * @return characterMap		HashMap<Character, String> 
	 */
	private HashMap<Character, String> setCharacterMapForRegularExpression(){

		HashMap<Character, String> characterMap = new HashMap<>();

		characterMap.put('A', "A");
		characterMap.put('C', "C");
		characterMap.put('G', "G");
		characterMap.put('T', "T");
		characterMap.put('R', "[AG]");
		characterMap.put('Y', "[CT]");
		characterMap.put('D', "[ATG]");
		characterMap.put('B', "[TGC]");
		characterMap.put('H', "[AUC]");
		characterMap.put('V', "[AGC]");
		characterMap.put('*', ".");
		return characterMap;
	}

	private HashMap<String, String> determineRefSeqIdsAssociatedToCoreProteinsAnnotatedByMotif(String annotationFile, String refSeqIdsByProteinFile){
		
		/* determine core proteins that are annotated by motif */
		HashSet<String> coreProteins = loadCoreProteinsForGivenMotif(annotationFile);

		/* get refSeqIds associated to proteins */
		HashMap<String, String> idMap = getRefSeqIds(refSeqIdsByProteinFile, coreProteins);

		return idMap;
	}

	private HashSet<String> loadCoreProteinsForGivenMotif(String annotationFile){

		HashSet<String> coreProteinSet = new HashSet<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			
			String line = in.readLine();
			while(line != null) {
				
				if(line.split("\t")[0].equals(this.motif)) {
					coreProteinSet.addAll(Arrays.asList(line.split("\t")[2].split("\\|")));
					break;
				}
				
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return coreProteinSet;
	}
	
	private HashMap<String, String> getRefSeqIds(String refSeqIdsByProteinFile, HashSet<String> coreProteins){
		
		HashMap<String, String> refSeqIdMap = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(refSeqIdsByProteinFile))));
			
			String line = in.readLine();
			
			while(line != null) {
				
				/* store refSeq IDs associated to core proteins annotated by motif  */
				String[] entry = line.split("\t");
				String prot = entry[0];
				
				/* check if protein is annotated by motif (ie. in set) */
				if(coreProteins.contains(prot)) {
					
					/* store info if protein contains sequences */
					if(entry.length > 1) {
						String[] ids = entry[1].split("\\|");	
						for(String id: ids) {
							refSeqIdMap.put(id, prot); 	// refSeq ID = protein name
						}
					}
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return refSeqIdMap;
	}

}
