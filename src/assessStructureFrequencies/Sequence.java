package assessStructureFrequencies;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import structures.Bulge;
import structures.End;
import structures.ExternalLoop;
import structures.Hairpin;
import structures.InnerLoop;
import structures.MultiLoop;
import structures.Stem;

public class Sequence {

	List<Stem> stems;
	List<Hairpin> hairpins;
	List<Bulge> bulges;
	List<InnerLoop> innerLoops;
	List<MultiLoop> multiLoops;
	List<ExternalLoop> externalLoops;
	List<End> ends;

	HashMap<String, Integer> frequencyMap;

	public Sequence(String file) {

		stems = new ArrayList<>();
		hairpins = new ArrayList<>();
		bulges = new ArrayList<>();
		innerLoops = new ArrayList<>();
		multiLoops = new ArrayList<>();
		externalLoops = new ArrayList<>();
		ends = new ArrayList<>();

		/* load structures */ 
		readFile(file);

		/* calculate frequencies */
		calculateFrequencies();

	}

	private void readFile(String file) {

		HashMap<Integer, Integer> multiLoopsMap = countMultiLoopSegments(file);

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(file))));

			String line = in.readLine();

			/* skip first 7 lines */ 
			for(int i=0; i<7; i++) {
				line = in.readLine();
			}

			
			int loopCount = 0;
			List<String> lines = new ArrayList<>();

			while(line != null) {

				/* load structures */ 
				switch(line.charAt(0)) {

				case 'S' : // Stem
					stems.add(new Stem(line));
					break;

				case 'H' : // Hairpin
					hairpins.add(new Hairpin(line));
					break;

				case 'B' : // Bulge
					bulges.add(new Bulge(line));
					break;

				case 'E': // End
					ends.add(new End(line));
					break;

				case 'X': // External loop
					externalLoops.add(new ExternalLoop(line));
					break;

				case 'I' : // Inner loop
					if(loopCount<2) {
						lines.add(line);
						loopCount++;
					}
					if(loopCount==2) {
						innerLoops.add(new InnerLoop(lines.get(0), lines.get(1)));
						lines.clear();
						loopCount=0;
					}
					break;

				case 'M': // Multi-loop
					String sub = line.substring(1, line.length());
					int current = Integer.parseInt(sub.split("\\.")[0]);

					if(loopCount < multiLoopsMap.get(current)) {
						lines.add(line);
						loopCount++;
					}

					if(loopCount == multiLoopsMap.get(current)) {
						multiLoops.add(new MultiLoop(lines));
						lines.clear();
						loopCount=0;
					}
					break;
				}
				line = in.readLine();

			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private HashMap<Integer, Integer> countMultiLoopSegments(String file){

		HashMap<Integer, Integer> multiLoopsMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(file))));

			String line = in.readLine();

			/* skip first 7 lines */ 
			for(int i=0; i<7; i++) {
				line = in.readLine();
			}

			while(line != null) {
				if(line.startsWith("M")) {

					String sub = line.substring(1, line.length());
					int current = Integer.parseInt(sub.split("\\.")[0]);

					if(multiLoopsMap.containsKey(current)) {
						multiLoopsMap.put(current, multiLoopsMap.get(current)+1);
					} else { 
						multiLoopsMap.put(current, 1);
					}
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return multiLoopsMap;
	}

	private void calculateFrequencies() {

		HashMap<String, Integer> freq = new HashMap<>();

		// stems
		for(Stem s : stems) {
			if(freq.containsKey(s.getStructureType())) {
				freq.put(s.getStructureType(), freq.get(s.getStructureType()) + 1);
			} else {
				freq.put(s.getStructureType(), 1);
			}
		}

		// hairpins
		for(Hairpin s : hairpins) {
			if(freq.containsKey(s.getStructureType())) {
				freq.put(s.getStructureType(), freq.get(s.getStructureType()) + 1);
			} else {
				freq.put(s.getStructureType(), 1);
			}
		}

		// bulges
		for(Bulge s : bulges) {
			if(freq.containsKey(s.getStructureType())) {
				freq.put(s.getStructureType(), freq.get(s.getStructureType()) + 1);
			} else {
				freq.put(s.getStructureType(), 1);
			}
		}

		// inner loop
		for(InnerLoop s : innerLoops) {
			if(freq.containsKey(s.getStructureType())) {
				freq.put(s.getStructureType(), freq.get(s.getStructureType()) + 1);
			} else {
				freq.put(s.getStructureType(), 1);
			}
		}

		// multi-loop
		for(MultiLoop s : multiLoops) {
			if(freq.containsKey(s.getStructureType())) {
				freq.put(s.getStructureType(), freq.get(s.getStructureType()) + 1);
			} else {
				freq.put(s.getStructureType(), 1);
			}
		}

		// external loops
		for(ExternalLoop s : externalLoops) {
			if(freq.containsKey(s.getStructureType())) {
				freq.put(s.getStructureType(), freq.get(s.getStructureType()) + 1);
			} else {
				freq.put(s.getStructureType(), 1);
			}
		}

		// end
		for(End s : ends) {
			if(freq.containsKey(s.getStructureType())) {
				freq.put(s.getStructureType(), freq.get(s.getStructureType()) + 1);
			} else {
				freq.put(s.getStructureType(), 1);
			}
		}
		
		this.frequencyMap = freq;
	}
	
	public List<Hairpin> getHairpins() {
		return hairpins;
	}
}
