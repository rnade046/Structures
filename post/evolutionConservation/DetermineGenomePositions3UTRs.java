package evolutionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;

public class DetermineGenomePositions3UTRs {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/Structures/post/nwTPD2/";
		String fastaFile = wd + "evolution/hg38_3utr_sequence.fasta.txt";
		String utrBedFile = wd + "evolution/evolutionConservation_genomicPositions_3UTR.bed";

		HashSet<String> utrPositions = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));			

			String line = in.readLine();

			while(line!=null) {

				/* check if sequence is associated to a protein in network */
				if(line.startsWith(">")) {

					//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
					String id = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];

					if(!line.contains("alt") && !line.contains("_fix")) {
						String coord = getUTRgenomicPosition(line);
						utrPositions.add(coord);
					}
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		printUTRBedFiles(utrPositions, utrBedFile);
	}


	private static String getUTRgenomicPosition(String header) {

		String chromosome = header.split(" ")[1].split("=")[1].split(":")[0];
		chromosome = chromosome+":";
		int startPos =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[0]);
		int endPos =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[1]);

		String genomicPosition = chromosome + startPos +"-"+endPos;

		return genomicPosition;
	}
	
	private static void printUTRBedFiles(HashSet<String> positions, String bedFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(bedFile)));

			out.write("Track - 3UTRs\n");

			for(String pos: positions) {

				String[] split = pos.split(":");
				out.write(split[0] +"\t"+ (Integer.parseInt(split[1].split("-")[0])) +"\t"+(Integer.parseInt(split[1].split("-")[1])) + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
