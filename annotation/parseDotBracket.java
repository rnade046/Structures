import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

public class parseDotBracket {

	public static void main(String[] args) {

		String inputFile = "/Users/rnadeau2/Documents/BEAR/test.txt";
		assessSecondaryStructure(inputFile);
	}

	public static void assessSecondaryStructure(String inputFile) {

		/* load sequence & dot bracket notation */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			char[] seq;
			char[] struct;
			int lineCount = 0;

			String line = in.readLine();
			while (line!=null) {



				if(lineCount==2) {
					struct = line.toCharArray();

					/* parse structure */
					parseStructure(struct);
					lineCount=0;
				}

				if(lineCount==1) {
					seq = line.toCharArray();
					lineCount++;
				}

				if(line.startsWith(">")) {
					lineCount++;
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void parseStructure(char[] structure) {

		/* search sequence for matching parenthesis */
		ArrayList<int[]> matchingBracketIdx = new ArrayList<>();

		for(int i=0; i<structure.length; i++) {

			if(structure[i] == '(') {
				int j = findClosingParen(structure, i);

				int[] currentIdx = new int[] {i, j};
				matchingBracketIdx.add(currentIdx);
			}
		}

		/* format structure */ 
		String[] formattedStructure = new String[structure.length];

		formattedStructure = findStems(matchingBracketIdx, formattedStructure);
		formattedStructure = findHairpins(matchingBracketIdx, structure, formattedStructure);
		formattedStructure = findMultiLoops(matchingBracketIdx, structure, formattedStructure);
		summarizeStems(formattedStructure);

		System.out.println("done");
	}

	private static int findClosingParen(char[] text, int openPos) {
		int closePos = openPos;
		int counter = 1;
		while (counter > 0) {
			char c = text[++closePos];
			if (c == '(') {
				counter++;
			}
			else if (c == ')') {
				counter--;
			}
		}
		return closePos;
	}

	/**
	 * Stems are paired base pairs; here we are looking to number the stems based on the identified matching brackets
	 * 
	 * @param matchingBracketIdx
	 * @param formattedStructure
	 * @return
	 */
	private static String[] findStems(ArrayList<int[]> matchingBracketIdx, String[] formattedStructure) {
		/* determine stems */
		int stemCount = 1;
		boolean newStem = true; 

		int currentStemLeft = 0;
		int currentStemRight = 0;

		int lastStemLeft = 1000;
		int lastStemRight = 0;

		for(int i=0; i<matchingBracketIdx.size(); i++) {

			currentStemLeft = matchingBracketIdx.get(i)[0];
			currentStemRight = matchingBracketIdx.get(i)[1];

			if(!newStem) { 

				/* check if same stem or different*/
				if(currentStemLeft == (lastStemLeft+1) && currentStemRight == (lastStemRight-1)) {
					formattedStructure[currentStemLeft] = "S" + stemCount;
					formattedStructure[currentStemRight] = "S" + stemCount;
				} else {
					newStem = true;
					stemCount++;
				}
			}

			if(newStem) {
				formattedStructure[currentStemLeft] = "S" + stemCount;
				formattedStructure[currentStemRight] = "S" + stemCount;

				newStem = false;
			}

			lastStemLeft = currentStemLeft;
			lastStemRight = currentStemRight;
		}
		return formattedStructure;
	}


	/**
	 * Hairpins are the single base pairs contained between brackets without interruption (...) 
	 * 
	 * @param matchingBracketIdx
	 * @param structure
	 * @param formattedStructure
	 * @return
	 */
	private static String[] findHairpins(ArrayList<int[]> matchingBracketIdx, char[] structure, String[] formattedStructure) {

		int hairpinCount = 1;
		for(int[] indexes: matchingBracketIdx) {

			boolean hairpin = true;
			for(int i=indexes[0]+1; i<indexes[1]; i++) {

				if(structure[i] != '.') {
					hairpin = false;
				}
			}

			if(hairpin) {
				for(int i=indexes[0]+1; i<indexes[1]; i++) {
					formattedStructure[i] = "H" + hairpinCount;
				}
				hairpinCount++;
			}
		}

		return formattedStructure;
	}

	private static String[] findMultiLoops(ArrayList<int[]> matchingBracketIdx, char[] structure, String[] formattedStructure) {

		List<int[]> mismatchBrackets = findMissMatchBrackets(structure);

		HashMap<Integer, Integer> matchingBrackets = findMatchingBrackets(mismatchBrackets, matchingBracketIdx);


		/* link miss match brackets */
		List<List<int[]>> linkedMissMatchedBrackets = new ArrayList<>();

		for(int i=0; i<mismatchBrackets.size(); i++) {
			int[] bracket = mismatchBrackets.get(i);

			if(bracket[0] + bracket[1] !=0) {

				List<int[]> currentLinks = new ArrayList<>();
				currentLinks.add(bracket);

				int openBracket = bracket[1];

				boolean linkBroken = false;
				while(!linkBroken) {

					if(matchingBrackets.containsKey(openBracket)) {

						int closingBracket = matchingBrackets.get(openBracket);
						boolean foundBracket = false;
						for(int j=i+1; j<mismatchBrackets.size(); j++) {

							int[] currentBracket = mismatchBrackets.get(j);
							if(currentBracket[0] == closingBracket) {
								currentLinks.add(currentBracket); // add to list
								openBracket = currentBracket[1];  // reset open bracket to search
								mismatchBrackets.set(j, new int[] {0, 0});
								foundBracket = true;
								break;
							}
						}

						if(!foundBracket) {
							currentLinks.add(new int[] {closingBracket, 0});
							linkedMissMatchedBrackets.add(currentLinks);
							linkBroken = true;
						}
					}
				}

			}
		}

		/* find final index */
		HashSet<Integer> loopsToRemove = new HashSet<>();
		for(int i=0; i<linkedMissMatchedBrackets.size(); i++) {

			List<int[]> loop = linkedMissMatchedBrackets.get(i);

			int currentIdx = loop.get(loop.size()-1)[0] +1;
			boolean indexFound = false;
			while(currentIdx < structure.length) {
				if(structure[currentIdx] == ')') {
					loop.get(loop.size()-1)[1] = currentIdx;
					indexFound = true;
					break;
				}
				currentIdx++;
			}

			if(!indexFound) {
				loopsToRemove.add(i);
			}
		}

		/* remove lists that aren't complete loops */
		List<List<int[]>> loops = new ArrayList<>();
		System.out.println("removed : ");
		for(int i=0; i<linkedMissMatchedBrackets.size(); i++) {
			if(!loopsToRemove.contains(i)) {
				loops.add(linkedMissMatchedBrackets.get(i));
			} else {
				System.out.print(i + ".");
			}
		}

		/* add initial two parenthesis; 0 = match of last, 1 = match of 2 (current 0) */
		/* determine indexes to match */
		HashMap<Integer, Integer> idxs = new HashMap<>();
		for(int i=0; i<loops.size(); i++) {

			List<int[]> loop = loops.get(i);
			idxs.put(loop.get(0)[0], 0);
			idxs.put(loop.get(loop.size()-1)[1], 0);
		}

		/* match indexes*/
		for(int[] matchingIdxs : matchingBracketIdx) {
			if(idxs.containsKey(matchingIdxs[0])) {
				idxs.replace(matchingIdxs[0], matchingIdxs[1]);
			}
			if(idxs.containsKey(matchingIdxs[1])) {
				idxs.replace(matchingIdxs[1], matchingIdxs[0]);
			}
		}

		/* set starting indexes */ 
		for(int i=0; i<loops.size(); i++) {

			List<int[]> loop = loops.get(i);

			int[] startOfLoop = new int[2];
			startOfLoop[0] = idxs.get(loop.get(loop.size()-1)[1]);
			startOfLoop[1] = idxs.get(loop.get(0)[0]);

			loops.get(i).add(0, startOfLoop);

		}

		/* set multi-loops */ 
		int loopCount = 1;
		for(List<int[]> l : loops) {
			for(int[] pos : l) {
				for(int i=pos[0]+1; i<pos[1]; i++) {
					if(structure[i] == '.') {
						formattedStructure[i] = "M" + loopCount;
					}
				}
			}
			loopCount++;
		}

		return formattedStructure;
	}

	private static void setLoopsAndBulges(ArrayList<int[]> matchingBracketIdx, char[] struct, String[] formattedStructure) {

		/* search for single strand */
		boolean ss = false; 
		int openBracket = 0;
		int closedBracket = 0;

		int countOpenNull = 0;
		int countClosedNull = 0;
		for(int i=0; i<formattedStructure.length; i++) {

			/* is single strand */
			if(formattedStructure[i] == null) {
				
				/* new single strand */
				if(!ss) {
					if(struct[i-1] == '(') {
						openBracket = i-1;
						ss = true;

						/* find closed bracket */
						for(int[] brackets : matchingBracketIdx) {
							if(brackets[0] == openBracket) {
								closedBracket = brackets[1];
								break;
							}
						}
						
						countOpenNull++;
					}
				}
				
			}
		}
	}

	private static List<int[]> findMissMatchBrackets(char[] structure){
		/* find parenthesis that are )...( */
		List<int[]> missMatchedParenthesis = new ArrayList<>();
		boolean potentialLoop = false;

		int openBracket = 0;
		int closedBracket = 0;

		for(int i=0; i<structure.length; i++) {

			if(potentialLoop) {
				if(structure[i] == '(') {
					openBracket = i;
					missMatchedParenthesis.add(new int[] {closedBracket, openBracket});
					potentialLoop = false;
				} 
				//				else if(structure[i] == ')') {
				//					potentialLoop = false;
				//				}
			} 

			if(structure[i] == ')') {
				closedBracket = i; 
				potentialLoop = true;
			}

		}
		return missMatchedParenthesis;
	}

	private static HashMap<Integer, Integer> findMatchingBrackets(List<int[]> mismatchBrackets, List<int[]> matchingBracketIdx){
		/* find closing indexes of open indexes for miss matched brackets */
		HashMap<Integer, Integer> matchingBrackets = new HashMap<>(); // open idx = closed idx for miss-matched brackets
		for(int[] missMatch : mismatchBrackets) {

			int openBracket = missMatch[1];
			for(int[] idxs : matchingBracketIdx) {
				if(idxs[0] == openBracket) {
					matchingBrackets.put(openBracket, idxs[1]);
					break;
				}
			}
		}
		return matchingBrackets;
	}


	private static void summarizeStems(String[] formattedStructure) {

		HashMap<String, Integer> stemCountMap = new HashMap<>();

		for(int i=0; i<formattedStructure.length; i++) {
			if(formattedStructure[i] != null) {

				if(stemCountMap.containsKey(formattedStructure[i])) {
					stemCountMap.replace(formattedStructure[i], stemCountMap.get(formattedStructure[i])+1);
				} else {
					stemCountMap.put(formattedStructure[i], 1);
				}
			}
		}

		HashMap<Integer, List<String>> stemMap = new HashMap<>();
		for(Entry<String, Integer> e : stemCountMap.entrySet()) {

			int baseNumbers = e.getValue()/2;

			if(stemMap.containsKey(baseNumbers)) {
				stemMap.get(baseNumbers).add(e.getKey());
			} else {
				List<String> listOfStems = new ArrayList<>();
				listOfStems.add(e.getKey());
				stemMap.put(baseNumbers, listOfStems);
			}
		}

		/* output stems */ 
		for(Entry<Integer, List<String>> e: stemMap.entrySet()) {
			System.out.println(e + " " + e.getValue().size());
		}
	}
}
