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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.Paranoid;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMReadGroupParser;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.ChromosomeSequence;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;

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
		final Paranoid paranoid = Paranoid.createThrowingInstance();
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
					if(tokens.length<21) throw new JvarkitException.TokenErrors(21, tokens);
					
					final SAMRecord rec = samRecordFactory.createSAMRecord(header);
					
					
					int col=0;
				    //0 matches - Number of bases that match that aren't repeats
				    final int matches = Integer.parseInt(tokens[col++]);
					//1 misMatches - Number of bases that don't match
				    final int misMatches = Integer.parseInt(tokens[col++]);
				    //2 repMatches - Number of bases that match but are part of repeats
				    final int repMatches=  Integer.parseInt(tokens[col++]);
				    //3 nCount - Number of "N" bases
				    final int nCount =  Integer.parseInt(tokens[col++]);
				    //4 qNumInsert - Number of inserts in query
				    final int qNumInsert =   Integer.parseInt(tokens[col++]);
				    //5 qBaseInsert - Number of bases inserted in query
				    final int qBaseInsert = Integer.parseInt(tokens[col++]);
				    //6 tNumInsert - Number of inserts in target
				    final int tNumInsert = Integer.parseInt(tokens[col++]);
				    //7 tBaseInsert - Number of bases inserted in target
				    final int tBaseInsert =  Integer.parseInt(tokens[col++]);
				    //8 strand - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
				    final char strand= tokens[col++].charAt(0);
				    //9 qName - Query sequence name
				    final String readName = tokens[col++];
				    //10 qSize - Query sequence size.
				    final int qSize = Integer.parseInt(tokens[col++]);
				    //11 qStart - Alignment start position in query
				    final int qStart =  Integer.parseInt(tokens[col++]);
				    //12 qEnd - Alignment end position in query
				    final int qEnd =  Integer.parseInt(tokens[col++]);
				    //13 tName - Target sequence name
				    final String tName = tokens[col++];
				    //14 tSize - Target sequence size
				    final int tSize =  Integer.parseInt(tokens[col++]);
				    //15 tStart - Alignment start position in target
				    final int tStart = Integer.parseInt(tokens[col++]);
				    //16 tEnd - Alignment end position in target
				    final int tEnd = Integer.parseInt(tokens[col++]);
				    //17 blockCount - Number of blocks in the alignment (a block contains no gaps)
				    final int blockCount = Integer.parseInt(tokens[col++]);
				    //18 blockSizes - Comma-separated list of sizes of each block. If the query is a protein and the target the genome, blockSizes are in amino acids. See below for more information on protein query PSLs.
					final int blockSizes[] = Arrays.stream(CharSplitter.COMMA.split(tokens[col++])).
							mapToInt(Integer::parseInt).
							toArray();
					//19
					final int qStarts[] = Arrays.stream(CharSplitter.COMMA.split(tokens[col++])).
							mapToInt(Integer::parseInt).
							toArray();
					//20
					final int tStarts[] = Arrays.stream(CharSplitter.COMMA.split(tokens[col++])).
							mapToInt(Integer::parseInt).
							toArray();

					final SAMSequenceRecord ssr = dict.getSequence(tName);
					if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(tName, dict);
					if(ssr.getSequenceLength()!=tSize) throw new IllegalArgumentException("expected contig leng="+ssr.getSequenceLength()+" for "+ssr.getSequenceName()+" but got "+tSize);
					
					rec.setReferenceName(tName);
					rec.setReadName(readName);
					rec.setAlignmentStart(tStart+1);
					rec.setMappingQuality(this.mapq);
					rec.setReadNegativeStrandFlag(strand=='-');
					rec.setAttribute("NM", misMatches);
					
					
					final List<CigarElement> cigar = new ArrayList<CigarElement>(blockCount+2);
					if(qStart>0) cigar.add(new CigarElement(qStart, CigarOperator.H));
					
					int y[]=qStarts;
					int z[]=tStarts;
					int gap_open =0;
					int gap_ext = 0;
					int readPos0 = y[0];
					int refPos0 = z[0];
					for(int i=1;i< blockCount;i++)
						{
						int q_stop = qStarts[i -1] + blockSizes[i-1];
						int t_stop = tStarts[i -1] + blockSizes[i-1];

						int q_gap_len = qStarts[i] - q_stop;
						int t_gap_len = tStarts[i] - t_stop;
						
						 ++gap_open;
						
						if (q_gap_len < t_gap_len) {
							  gap_ext += t_gap_len - q_gap_len;
							  cigar.add(new CigarElement(y[i] - readPos0 ,CigarOperator.M));
							  cigar.add(new CigarElement(t_gap_len - q_gap_len,CigarOperator.D));
							  readPos0 = qStarts[i];
							  refPos0 = tStarts[i];
							} else if (t_gap_len < q_gap_len) { //
							  gap_ext += q_gap_len - t_gap_len;
							  cigar.add(new CigarElement(z[i] - refPos0,CigarOperator.M));
							  cigar.add(new CigarElement(q_gap_len - t_gap_len,CigarOperator.I));
							  readPos0 = qStarts[i];
							  refPos0 = tStarts[i];
							} else
							{
								//nothing
							}
						}
					cigar.add(new CigarElement(qEnd - readPos0,CigarOperator.M));

					if (qSize != qEnd) cigar.add(new CigarElement(qSize - qEnd ,CigarOperator.H));
					
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
						final StringBuilder readseq = new StringBuilder();
						int ref0= tStart;
						for(final CigarElement ce: cigar) {
							switch(ce.getOperator()) {
								case H:break;
								case D: case N:ref0+=ce.getLength();break;
								case I: readseq.append(StringUtils.repeat(ce.getLength(),'N'));break;
								case M: case EQ: case X:
									{
									readseq.append(genomicSeq.subSequence(ref0, ref0+ce.getLength())).toString().toUpperCase();
									ref0+=ce.getLength();	
									break;
									}
								default: throw new IllegalStateException(""+ce.getOperator());
								}
							
							}
						final String qual = StringUtils.repeat(readseq.length(), this.basequality);  
						rec.setReadString(readseq.toString());
						rec.setBaseQualityString(qual);
						}
					if(readGroupId!=null) rec.setAttribute("RG",readGroupId);
					rec.setCigar(new Cigar(cigar));
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
