/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMReadGroupParser;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.ChromosomeSequence;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.LineIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.samtools.GenomeSequence;

/**

BEGIN_DOC

## Motivation

Convert **OPS** to **PSL** http://genome.ucsc.edu/FAQ/FAQformat.html#format2 or **BED12** .


END_DOC
*/
@Program(name="pslx2bam",
	deprecatedMsg="use bedtools/bamtobed",
	description="Convert PSLX to SAM/BAM",
	keywords={"sam","bam","psl","pslx"}
	)
public class PslxToBam extends Launcher
	{
	private static final Logger LOG = Logger.build(PslxToBam.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"-g","--rg"},description= SAMReadGroupParser.OPT_DESC)
	private String readGroupStr = null;
	@Parameter(names= {"-q","--mapq"},description="default mapping quality")
	private int mapq = SAMRecord.UNKNOWN_MAPPING_QUALITY;
	@Parameter(names= {"-b","--baq"},description="default base symbol")
	private char basequality = '#';
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	
	@Override
	public int doWork(final List<String> args)
		{
		SAMFileWriter sw=null;
		BufferedReader br=null;
		ReferenceSequenceFile referenceSequenceFile=null;
		ChromosomeSequence genomicSeq=null;
		try
			{
			referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx);
			
			br = super.openBufferedReader(oneFileOrNull(args));
			this.writingBamArgs.setReferencePath(this.faidx);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(referenceSequenceFile);
			final SAMFileHeader header= new SAMFileHeader(dict);
			header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
			final String readGroupId;
			if(!StringUtils.isBlank(readGroupStr)) {
				final SAMReadGroupParser rgParser = new SAMReadGroupParser();
				final SAMReadGroupRecord rg =rgParser.apply(this.readGroupStr);
				header.addReadGroup(rg);
				readGroupId = rg.getId();
				}
			else
				{
				readGroupId=null;
				}
			JVarkitVersion.getInstance().addMetaData(this, header);
			final SAMRecordFactory samRecordFactory = DefaultSAMRecordFactory.getInstance();
			
			sw = writingBamArgs.openSamWriter(this.outputFile, header, true);
			
			String line;
			while((line=br.readLine())!=null)
					{
					if(StringUtils.isBlank(line) || 
		    				line.startsWith("----") ||
		    				line.startsWith("psLayout") || 
		    				!Character.isDigit(line.charAt(0)))//PSL header, empty lines...
		    				{
		    				continue;
		    				}
					final String tokens[] = CharSplitter.TAB.split(line);
					
					
					final SAMRecord rec = samRecordFactory.createSAMRecord(header);
					
					
					String readName = tokens[9];
					
					int match= Integer.parseInt(tokens[0]);
					int mismatch= Integer.parseInt(tokens[1]);
					int repmatch= Integer.parseInt(tokens[2]);
					int countN= Integer.parseInt(tokens[3]);
					
					int querySize = Integer.parseInt(tokens[12]);
					int queryStart = Integer.parseInt(tokens[13]);
					int queryEnd = Integer.parseInt(tokens[14]);
					
					boolean negativeStrand = tokens[8].equals("-");
					rec.setReferenceName(tokens[13]);
					rec.setReadName(readName);
					if(negativeStrand) {
						rec.setAlignmentStart(Integer.parseInt(tokens[15])+1);
						}
					else
						{
						rec.setAlignmentStart(Integer.parseInt(tokens[16])+1);
						}
					rec.setMappingQuality(this.mapq);
					rec.setReadNegativeStrandFlag(negativeStrand);
					rec.setAttribute("NM", mismatch);
					
					final int blockCount = Integer.parseInt(tokens[17]);
					
					int blockSizes[] = Arrays.stream(CharSplitter.COMMA.split(tokens[18])).
							mapToInt(Integer::parseInt).
							toArray();
					int qStarts[] = Arrays.stream(CharSplitter.COMMA.split(tokens[19])).
							mapToInt(Integer::parseInt).
							toArray();
					int tStarts[] = Arrays.stream(CharSplitter.COMMA.split(tokens[20])).
							mapToInt(Integer::parseInt).
							toArray();
					final StringBuilder cigar = new StringBuilder();
					int gap_open =0;
					int gap_ext = 0;
					int y0 = qStarts[0];
					int z0 = tStarts[0];
					for(int i=1;i< blockCount;i++)
						{
						int ly = qStarts[i] - qStarts[i -1] - blockSizes[i-1];
						int lz = tStarts[i] - tStarts[i-1] -blockSizes[i-1];
						if (ly < lz) {
							  ++gap_open;
							  gap_ext += lz - ly;
							  cigar.append(qStarts[i] - y0).append("M");
							  cigar.append(lz - ly).append("D");
							  y0 = qStarts[i];
							  z0 = tStarts[i];
							} else if (lz < ly) { //
							  ++gap_open;
							  gap_ext +=ly - lz;
							  cigar.append(tStarts[i] - z0).append("M");
							  cigar.append(lz - lz).append("D");
							  y0 = qStarts[i];
							  z0 = tStarts[i];
							}
						
						}
					
					cigar.append(0).append("M");
					if (tokens[10] != tokens[12]) cigar.append(0).append('H') ;
					
					if(tokens.length>22) {
						final String readseq = tokens[22].replace(",", "").trim().toUpperCase();
						final String qual = StringUtils.repeat(readseq.length(), this.basequality);  
						rec.setReadString(readseq);
						rec.setBaseQualityString(qual);
						}
					else
						{
						if(genomicSeq==null || !genomicSeq.getChrom().equals(tokens[13])) {
							genomicSeq = new GenomicSequence(referenceSequenceFile,tokens[13]);
							}
						//TODO
						}
					if(readGroupId!=null) rec.setAttribute("RG",readGroupId);
					rec.setCigarString(cigar.toString());
					sw.addAlignment(rec);
					}
				
			sw.close();sw=null;
			br.close();br=null;
			referenceSequenceFile.close();
			referenceSequenceFile=null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(referenceSequenceFile);
			CloserUtil.close(sw);
			CloserUtil.close(br);
			}
		}
	public static void main(final String[] args)
		{
		new PslxToBam().instanceMainWithExit(args);
		}
	}
