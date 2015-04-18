import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.*;
import java.util.HashMap;
import java.util.regex.*;


public class SNVSupportReads {
	public static void main(String[] args) throws IOException {
		String bamFilePath = "test/normal.chr11_533779_533779.bam";
		int interestedPosition = 533779 ;
		char alt = 'C';
		new SNVSupportReads().BaseInPos(new File(bamFilePath), interestedPosition, alt);
		
		//String bamFilePath = args[0];
		//int interestedPosition = Integer.parseInt(args[1]);
		//char alt = args[2].toCharArray()[0];
		//new SNVSupportReads().BaseInPos(new File(bamFilePath), interestedPosition, alt);
	}
	
	public void BaseInPos(File bamFile, int interestedPos, char alt) throws IOException {
		int totalValidateReads = 0;
		int totalSNVSupportValidateReads = 0;
		
		final SamReader inputSam = SamReaderFactory.makeDefault().open(bamFile); 
		for (final SAMRecord samRecord : inputSam) {
			
			//if(samRecord.getMappingQuality() == 0) {
			//	continue;
			//}
			
			if(interestedPos > samRecord.getAlignmentEnd() || interestedPos < samRecord.getAlignmentStart()) {
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
					if(posRef == interestedPos && samRecord.getBaseQualities()[posSEQ] >= 13) {
						if(samRecord.getReadString().toCharArray()[posSEQ] == alt) {
							totalSNVSupportValidateReads++;
						}
						totalValidateReads++;
						//System.out.println(samRecord.getReadString().toCharArray()[posSEQ]);
					}
					posRef++;
					posSEQ++;
				} else if(op == 'I') {
					posSEQ++;
				} else if(op == 'D') {
					posRef++;
				}

				idx_nonClipCigarExtArr++;
			}
		}
		inputSam.close();
		
		// TMP
		//System.out.println(totalValidateReads + "\t" + totalSNVSupportValidateReads + "\t" + (double)totalSNVSupportValidateReads/(double)totalValidateReads);
		System.out.print(totalValidateReads + "\t" + totalSNVSupportValidateReads + "\t" + (double)totalSNVSupportValidateReads/(double)totalValidateReads + "\t");
	}
}
