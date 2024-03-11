package investigateModules;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class Annotation {

	private List<String> proteins;
	private List<Double> bpScores;
	private List<Double> percentiles;

	public Annotation(String[] annotationProteinInfo, String[] percentileProteinInfo) {
		
		this.proteins = getAnnotatedProteins(annotationProteinInfo);
		this.bpScores = getScores(annotationProteinInfo);
		this.percentiles = getScores(percentileProteinInfo);
	}
	
	private List<String> getAnnotatedProteins(String[] values){
		
		List<String> annotatedProteins = new ArrayList<>();
		
		for(String v : values) {
			annotatedProteins.add(v.split("\\_")[0]);
		}
		return annotatedProteins;
 	}

	private List<Double> getScores(String[] values){
		
		List<Double> scores = new ArrayList<>();
		
		for(String v : values) {
			scores.add(Double.parseDouble(v.split("\\_")[1]));
		}
		return scores;
	}

	public List<String> getProteinList(){
		return this.proteins;
	}
	
	public List<Double> getPercentileScores(){
		return this.percentiles;
	}

	public void assessDegreeDistribution(HashMap<String, Protein> proteinMap, String outputFile) {

		/* get list of degrees for subset of proteins */		
		List<Integer> degrees = new ArrayList<>();

		for(String prot: this.proteins) {
			degrees.add(proteinMap.get(prot).getDegree());
		}

		/* count frequency of degrees */
		int max = Collections.max(degrees);

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("degree\tcount\tnormalizeCount\n");
			for(int i=0; i<=max; i=i+10) {

				int countDegrees = 0;
				for(int d: degrees) {
					if(d >= i && d < (i+10)) {
						countDegrees++;
					}
				}
				out.write(i + "\t" + countDegrees + "\t" + (countDegrees / (double) this.proteins.size()) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void assessBPscoreDistribution(String outputFile) {
		
		double min = Math.floor(Collections.min(this.bpScores));
		double max = Math.ceil(Collections.max(this.bpScores));
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("score\tcount\tnormalizeCount\n");
			for(double i=min; i<max; i++) {

				int countScores = 0;
				for(double s: this.bpScores) {
					if(s >= i && s < (i+1)) {
						countScores++;
					}
				}
				out.write(i + "\t" + countScores + "\t" + (countScores / (double) this.proteins.size()) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

public void assessPercentileScoreDistribution(String outputFile) {
		
		double min = Math.floor(Collections.min(this.percentiles));
		double max = Math.ceil(Collections.max(this.percentiles));
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("percentile\tcount\tnormalizeCount\n");
			for(double i=min; i<max; i++) {

				int countScores = 0;
				for(double s: this.percentiles) {
					if(s >= i && s < (i+1)) {
						countScores++;
					}
				}
				out.write(i + "\t" + countScores + "\t" + (countScores / (double) this.proteins.size()) + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
