import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class INDELSupportReads {

	public static void main(String[] args) throws IOException {
		
		// TMP
		//BufferedReader brNLH = new BufferedReader(new FileReader("test/patient_nrt_low_high"));
		//String lineNLH = null;
		//HashMap<String, String[]> map = new HashMap<String, String[]>();
		//while((lineNLH = brNLH.readLine()) != null) {
		//	String[] arr = lineNLH.split("\t");
		//	map.put(arr[0], new String[]{arr[1], arr[2], arr[3]});
		//}
		//brNLH.close();
		
		//BufferedReader br = new BufferedReader(new FileReader("test/INDEL_coords.txt"));
		//String line = null;
		//while((line = br.readLine()) != null) {
		//	String[] arr = line.split("\t");
		//	String p = arr[0];
		//	String chr = arr[1];
		//	int pos = Integer.parseInt(arr[2]);
		//	String ref = arr[4];
		//	String alt = arr[5];
		//	
		//	char type = ' ';
		//	String indelString = "";
		//	if(ref.equals("-")) {
		//		type = 'I';
		//		indelString = alt;
		//	} else {
		//		type = 'D';
		//		indelString = ref;
		//	}
		//	
		//	new INDELSupportReads().IndelInPos(new File("test/" + map.get(p)[2] + "_" + chr + ".bam"), pos, type, indelString); // HIGH_WGA
		//	new INDELSupportReads().IndelInPos(new File("test/" + map.get(p)[1] + "_" + chr + ".bam"), pos, type, indelString); // LOW_WGA
		//	new INDELSupportReads().IndelInPos(new File("test/" + map.get(p)[0] + "_" + chr + ".bam"), pos, type, indelString); // NRT
		//	System.out.println();
		//}
		//br.close();
		
		String bamFilePath = args[0];
		int interestedPosition = Integer.parseInt(args[1]);
		char type = args[2].toCharArray()[0];
		String indelString = args[3];
		new INDELSupportReads().IndelInPos(new File(bamFilePath), interestedPosition, type, indelString);
	}
	
	public void IndelInPos(File bamFile, int interestedPos, char type/*I or D*/, String indelString) throws IOException {
		
		int totalValidateReads = 0;
		int totalIndelSupportValidateReads = 0;
		
		if( type == 'I') {
			interestedPos += 1;
		}
		
		final SamReader inputSam = SamReaderFactory.makeDefault().open(bamFile); 
		for (final SAMRecord samRecord : inputSam) {
			
			if(interestedPos < samRecord.getAlignmentStart() || interestedPos > samRecord.getAlignmentEnd()) {
				continue;
			}
			
			totalValidateReads++;
		
			if(samRecord.getMappingQuality() == 0) {
				continue;
			}
			
			String cigarString = samRecord.getCigarString();
			Pattern pHeadSP = Pattern.compile("^(\\d+H)*((\\d+)S)+.*");
			Matcher mHeadSP = pHeadSP.matcher(cigarString);
			int posSEQ = 0;
			if(mHeadSP.matches()) {
				posSEQ = Integer.parseInt(mHeadSP.group(3)); // not count soft clip bases
			}
			
			Pattern nonClipCigarPattern = Pattern.compile("(\\d+)([MID])");
			Matcher nonClipCigarMatcher = nonClipCigarPattern.matcher(cigarString);
			String nonClipCigarExt = "";
			while(nonClipCigarMatcher.find()) {
				int n = Integer.parseInt(nonClipCigarMatcher.group(1));
				String op = nonClipCigarMatcher.group(2);
				nonClipCigarExt += new String(new char[n]).replace("\0", op);
			}
			
			char[] nonClipCigarExtArr = nonClipCigarExt.toCharArray();
			int posRef = samRecord.getAlignmentStart();
			int idx_nonClipCigarExtArr = 0;
			
			while(posRef <= interestedPos && idx_nonClipCigarExtArr < nonClipCigarExtArr.length) {
				char op = nonClipCigarExtArr[idx_nonClipCigarExtArr];
				if(op == 'M') {
					posRef++;
					posSEQ++;
				} else if(op == 'I') {
					if(posRef == interestedPos) {
						int indelLen = indelString.length();
						String indelCigar = new String(Arrays.copyOfRange(nonClipCigarExtArr, idx_nonClipCigarExtArr, idx_nonClipCigarExtArr + indelLen));
						String imagineCigar = new String(new char[indelLen]).replace("\0", "I");
						if (indelCigar.equals(imagineCigar)) {
							String insertionSEQ = samRecord.getReadString().substring(posSEQ, posSEQ + indelLen);
							if(insertionSEQ.equals(indelString)) {
								totalIndelSupportValidateReads++;
								//System.out.println(samRecord.getReadName() + "\t" + samRecord.getReferenceName() + "\t" + samRecord.getAlignmentStart() + "\t" +samRecord.getCigarString() + "\t" + samRecord.getReadString());
							}
						}
						break;
					}
					posSEQ++;
				} else if(op == 'D') {
					if(posRef == interestedPos) {
						int indelLen = indelString.length();
						String indelCigar = new String(Arrays.copyOfRange(nonClipCigarExtArr, idx_nonClipCigarExtArr, idx_nonClipCigarExtArr + indelLen));
						String imagineCigar = new String(new char[indelLen]).replace("\0", "D");
						if (indelCigar.equals(imagineCigar)) {
							totalIndelSupportValidateReads++;
							//System.out.println(samRecord.getReadName() + "\t" + samRecord.getReferenceName() + "\t" + samRecord.getAlignmentStart() + "\t" +samRecord.getCigarString() + "\t" + samRecord.getReadString());
						}
						break;
					}
					posRef++;
				}
				idx_nonClipCigarExtArr++;
			}
		}
		inputSam.close();
		
		// TMP
		//System.out.println(totalValidateReads + "\t" + totalIndelSupportValidateReads + "\t" + (double)totalIndelSupportValidateReads/(double)totalValidateReads);
		System.out.print(totalValidateReads + "\t" + totalIndelSupportValidateReads + "\t" + (double)totalIndelSupportValidateReads/(double)totalValidateReads + "\t");
	}
}