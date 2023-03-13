package structures;

import java.util.ArrayList;
import java.util.List;

import assessStructureFrequencies.Sequence;

public class Segment {

	String id;
	int[][] bounds;
	List<String> seqs;
	int length;
	String structureType;

	List<Stem> stems;
	Hairpin hairpin;
	List<Bulge> bulges;
	List<InnerLoop> innerLoops;
	List<MultiLoop> multiLoops;
	List<ExternalLoop> externalLoops;
	List<End> ends;

	public Segment(String line, String seq, Sequence s) {
		
		stems = new ArrayList<>();
		bulges = new ArrayList<>();
		innerLoops = new ArrayList<>();
		multiLoops = new ArrayList<>();
		externalLoops = new ArrayList<>();
		ends = new ArrayList<>();

		String[] elements = line.split("\\s+");
		id = elements[0];

		bounds = new int[2][2];
		String[] start = elements[2].split("\\..");
		String[] end =  elements[4].split("\\..");

		bounds[0] = new int[] {Integer.parseInt(start[0]), Integer.parseInt(start[1])};
		bounds[1] = new int[] {Integer.parseInt(end[0]), Integer.parseInt(end[1])};
		
		length = Integer.parseInt(elements[1].split("bp")[0]);

		
	}
	
	private void searchForStructures(String seq, Sequence s) {
		
		/* get substrings */
		char[] s1 = seq.substring(bounds[0][0], bounds[0][1]).toCharArray();
		boolean[] s1Found = new boolean[s1.length];
		
		char[] s2 = seq.substring(bounds[1][0], bounds[1][1]).toCharArray();
		boolean[] s2Found = new boolean[s2.length];
		
		/* hairpin index */
		char lower = seq.charAt(bounds[0][1] + 1);
		char upper = seq.charAt(bounds[1][0] - 1);
		
		if(lower == 'H' && upper == 'H') {
			
			int[] hairpinBounds = new int[] {bounds[0][1] + 1, bounds[1][0] -1};
			hairpin = findHairpin(hairpinBounds, s.getHairpins());
		}
		
		
		
	}
	
	private Hairpin findHairpin(int[] hairpinBounds, List<Hairpin> hairpinList) {
		
		for(Hairpin h : hairpinList) {
			if(hairpinBounds[0] == h.bounds[0]) {
				return h;
			}
		}
	}
	 
}
