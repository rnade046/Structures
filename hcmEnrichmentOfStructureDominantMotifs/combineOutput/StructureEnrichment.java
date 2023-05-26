package combineOutput;

import java.util.List;

public class StructureEnrichment {

	private char structureType;
	private List<Double> significantScore;
	private List<Double> enrichmentScore;

	public StructureEnrichment(char s, List<Double> sScores, List<Double> eScores) {
		this.structureType = s;
		this.significantScore = sScores;
		this.enrichmentScore = eScores;
	}
	
	public char getStructureType() {
		return this.structureType;
	}
	
	public double[] getFormattedEnrichmentScores() {
		
		double[] formattedEnrichmentScores = new double[enrichmentScore.size()];
		
		for(int i=0; i<significantScore.size(); i++) {
			
			if(significantScore.get(i) <= 0.005) {
				formattedEnrichmentScores[i] = enrichmentScore.get(i);
			}
			
		}
		return formattedEnrichmentScores;
	}
}
