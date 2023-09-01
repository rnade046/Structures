package sampling;

import java.util.ArrayList;
import java.util.List;

public class Protein {
	
	private String proteinName;
	private int annotationFrequncy;
	private List<Double> scores;
	
	public Protein(String name, double score) {
		this.proteinName = name;
		this.annotationFrequncy = 1;
		
		this.scores = new ArrayList<>();
		this.scores.add(score);
	}
	
	public String getProteinName() {
		return this.proteinName;
	}
	
	public int getAnnotationFrequency() {
		return this.annotationFrequncy;
	}
	
	public List<Double> getScores(){
		return this.scores;
	}
	
	public void updateFrequency() {
		this.annotationFrequncy++;
	}
	
	public void addScores(double s) {
		this.scores.add(s);
	}
}
