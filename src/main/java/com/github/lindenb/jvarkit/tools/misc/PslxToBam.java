/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
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
import htsjdk.samtools.SAMProgramRecord;
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

## motivation

convert the psl  (or pslx) format to sam/bam.

When the input format is pslx. The program will use the bases to set the SEQ column. Nevertheless, there are sometimes fewers bases in the pslx format than expected (I don't understand why), so I sometimes fill the SEQ with 'N'.

## Example

### remapping reads with blat

```
samtools fasta src/test/resources/S1.bam |\
	blat -out=pslx  src/test/resources/rotavirus_rf.fa stdin stdout |\
	java -jar dist/psl2bam.jar -R src/test/resources/rotavirus_rf.fa 
```

### example
 
```
java -jar dist/psl2bam.jar -R src/test/resources/rotavirus_rf.fa   input.psl

@HD	VN:1.6	SO:unsorted
@SQ	SN:RF01	LN:3302	M5:59dccb944425dd61f895a564ad7b56a7	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
(...)
@SQ	SN:RF11	LN:666	M5:7a7cf2c7813f2e8bd74be383014202ca	UR:https://raw.githubusercontent.com/lindenb/jvarkit/master/src/test/resources/rotavirus_rf.fa	SP:rotavirus
@PG	ID:0	CL:-R src/test/resources/rotavirus_rf.fa	PN:psl2bam	VN:a08d9c1
@CO	psl2bam. compilation:20190918074721 githash:a08d9c1 htsjdk:2.20.1 date:20190918075100. cmd:-R src/test/resources/rotavirus_rf.fa
RF01:100-200	0	RF01	100	255	101M	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:0
RF01:100-200/rc	16	RF01	100	255	101M	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:0
RF01:100-200+N	0	RF01	100	255	5H101M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:0
RF01:100-200/rc+N	16	RF01	100	255	5H101M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:0
RF01:100-200+N+mismatch	0	RF01	100	255	5H101M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:1
RF01:100-200/rc+N+mismatch	16	RF01	100	255	5H101M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:1
RF01:100-200+N+DEL	0	RF01	100	255	5H33M29D39M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTCTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:29
RF01:100-200/rc+N+DEL	16	RF01	100	255	5H33M29D39M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTCTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:29
RF01:100-200+N+INS	0	RF01	100	255	5H31M29I70M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAAAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:29
RF01:100-200/rc+N+INS	16	RF01	100	255	5H67M29I34M5H	*	0	0	TATTCTTCCAATAGTGAATTAGAGAATAGATGTATTGAATTTCATTCTAAATGCTTAGAAAACTCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGAATGGACTATCATTGAAAAAGCTCTTTGTTGA	2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222	PG:Z:0	NM:i:29
```

## see also

 * https://github.com/samtools/samtools/blob/develop/misc/psl2sam.pl
 * https://github.com/bsipos/uncle_psl

END_DOC
*/
@Program(name="psl2bam",
	description="Convert PSL to SAM/BAM",
	keywords={"blat","sam","bam","psl","pslx"},
	creationDate="20190918",
	modificationDate="20190918"
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
	private char basequality = '2';
	@Parameter(names= {"-D","--disable-secondary"},description="disable 'consecutive reads with same name will be flagged as secondary'")
	private boolean disable_secondary_flag=false;
	@Parameter(names= {"--disable-pslx"},description="disable use of bases provided by the pslX format.")
	private boolean disable_pslx_bases=false;
	@Parameter(names= {"--insert-base"},description="default base for insertion")
	private char insert_base='N';
	@Parameter(names= {"--failed"},description="save problematic lines here.",hidden=true)
	private Path failedPath = null;
	@Parameter(names= {"--intron"},description="use 'N' operator instead of 'D' when deletion are larger than 'x'")
	private int use_N_operator_size = 50;

	
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
			final SAMProgramRecord prg = header.createProgramRecord();
			prg.setCommandLine(super.getProgramCommandLine());
			prg.setProgramName(getProgramName());
			prg.setProgramVersion(getGitHash());
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			final SAMRecordFactory samRecordFactory = DefaultSAMRecordFactory.getInstance();
			
			sw = writingBamArgs.openSamWriter(this.outputFile, header, true);
			
			
			PrintWriter failed;
			if(this.failedPath==null) {
				failed = new PrintWriter(new NullOuputStream());
			} else
				{
				failed = IOUtils.openPathForPrintWriter(failedPath);
				}
			
			
			String prevReadName=null;
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
					if(tokens.length<21) {
						sw.close();
						br.close();
						throw new JvarkitException.TokenErrors(21, tokens);
					}
					
					final SAMRecord rec = samRecordFactory.createSAMRecord(header);
					
					
					int col=0;
				    //0 matches - Number of bases that match that aren't repeats
				    @SuppressWarnings("unused")
					final int matches = Integer.parseInt(tokens[col++]);
					//1 misMatches - Number of bases that don't match
				    final int misMatches = Integer.parseInt(tokens[col++]);
				    //2 repMatches - Number of bases that match but are part of repeats
				    @SuppressWarnings("unused")
				    final int repMatches=  Integer.parseInt(tokens[col++]);
				    //3 nCount - Number of "N" bases
				    @SuppressWarnings("unused")
				    final int nCount =  Integer.parseInt(tokens[col++]);
				    //4 qNumInsert - Number of inserts in query
				    @SuppressWarnings("unused")
				    final int qNumInsert =   Integer.parseInt(tokens[col++]);
				    //5 qBaseInsert - Number of bases inserted in query
				    @SuppressWarnings("unused")
				    final int qBaseInsert = Integer.parseInt(tokens[col++]);
				    //6 tNumInsert - Number of inserts in target
				    @SuppressWarnings("unused")
				    final int tNumInsert = Integer.parseInt(tokens[col++]);
				    //7 tBaseInsert - Number of bases inserted in target
				    @SuppressWarnings("unused")
				    final int tBaseInsert =  Integer.parseInt(tokens[col++]);
				    //8 strand - "+" or "-" for query strand. For translated alignments, second "+"or "-" is for target genomic strand.
				    final char strand= tokens[col++].charAt(0);
				    //9 qName - Query sequence name
				    final String readName = tokens[col++];
				    //10 qSize - Query sequence size.
				    final int qSize = Integer.parseInt(tokens[col++]);
				    //11 qStart - Alignment start position in query
				    int qStart =  Integer.parseInt(tokens[col++]);
				    //12 qEnd - Alignment end position in query
				    int qEnd =  Integer.parseInt(tokens[col++]);
				    //13 tName - Target sequence name
				    final String tName = tokens[col++];
				    //14 tSize - Target sequence size
				    final int tSize =  Integer.parseInt(tokens[col++]);
				    //15 tStart - Alignment start position in target
				    final int tStart = Integer.parseInt(tokens[col++]);
				    //16 tEnd - Alignment end position in target
				    @SuppressWarnings("unused")
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

					if (strand=='-') {
					    int tmp = qStart;
					    qStart = qSize - qEnd;
					    qEnd = qSize - tmp;
					  }

					
					
					final SAMSequenceRecord ssr = dict.getSequence(tName);
					if(ssr==null) {
						throw new JvarkitException.ContigNotFoundInDictionary(tName, dict);
					}
					if(ssr.getSequenceLength()!=tSize) throw new IllegalArgumentException("expected contig leng="+ssr.getSequenceLength()+" for "+ssr.getSequenceName()+" but got "+tSize);
					
					rec.setReferenceName(tName);
					rec.setReadName(readName);
					rec.setAlignmentStart(tStart+1);
					rec.setMappingQuality(this.mapq);
					rec.setReadNegativeStrandFlag(strand=='-');
					
					
					final List<CigarElement> cigar = new ArrayList<CigarElement>(blockCount*2+2);
					if(qStart>0) cigar.add(new CigarElement(qStart, CigarOperator.H));
					
					int gap_ext = 0;
					
					try {
						int readPos0 = qStarts[0];
						int refPos0 = tStarts[0];
						for(int i=1;i< blockCount;i++)
							{
							final int q_stop = qStarts[i -1] + blockSizes[i-1];
							final int t_stop = tStarts[i -1] + blockSizes[i-1];
	
							final int q_gap_len = qStarts[i] - q_stop;
							final int t_gap_len = tStarts[i] - t_stop;						
							
							if (q_gap_len < t_gap_len) {
								 final int cigar_size = t_gap_len - q_gap_len;
								  gap_ext += cigar_size;
								  cigar.add(new CigarElement(qStarts[i] - readPos0 ,CigarOperator.M));
								  cigar.add(new CigarElement(cigar_size,cigar_size> this.use_N_operator_size ?CigarOperator.N:CigarOperator.D));
								  readPos0 = qStarts[i];
								  refPos0 = tStarts[i];
								} 
							else if (t_gap_len < q_gap_len) {
								  final int cigar_size = q_gap_len - t_gap_len;
								  gap_ext += cigar_size;
								  cigar.add(new CigarElement(tStarts[i] - refPos0,CigarOperator.M));
								  cigar.add(new CigarElement(cigar_size,CigarOperator.I));
								  readPos0 = qStarts[i];
								  refPos0 = tStarts[i];
								} 
							else
								{
									//nothing
								}
							}
						cigar.add(new CigarElement(qEnd - readPos0,CigarOperator.M));
	
						if (qSize > qEnd) cigar.add(new CigarElement(qSize - qEnd ,CigarOperator.H));
						
					} catch(final IllegalArgumentException err) {
						LOG.error("Cannot build cigar string for \""+line+"\"",err);
						failed.append(line);
						continue;
						}
					
					
					if(!disable_secondary_flag && readName.equals(prevReadName)) {
						rec.setSecondaryAlignment(true);
					}
					prevReadName = readName;

					rec.setAttribute("NM", misMatches);
					rec.setAttribute("AS", qSize-(misMatches + gap_ext));

					
					final String readBases;
					final int readbases_column_index=21;
					// tokens[21] is QUERY sequence in nucleic acid when using pslx format
					if(!disable_pslx_bases && tokens.length> readbases_column_index ) {
						final StringBuilder readBaseBuilder = new StringBuilder(tokens[readbases_column_index].length());
						final String t2[]  = CharSplitter.COMMA.split(tokens[readbases_column_index]);
						if(t2.length!=blockCount) throw new JvarkitException.TokenErrors(blockCount,t2);
						for(int j=0;j< t2.length;++j) {
							if(t2[j].length()!=blockSizes[j])  throw new IllegalAccessException("t2["+j+"]="+t2[j]+" expected length;"+blockSizes[j]);
							readBaseBuilder.append(t2[j].toUpperCase());
							}
						readBases = readBaseBuilder.toString();
						}
					else
						{
						readBases=null;
						}
					
						{
						if(genomicSeq==null || !genomicSeq.getChrom().equals(tName)) {
							genomicSeq = new GenomicSequence(referenceSequenceFile,tName);
							}
						final StringBuilder readseq = new StringBuilder();
						int ref0= tStart;
						int read0=0;
						for(final CigarElement ce: cigar) {
							switch(ce.getOperator()) {
								case H:break;
								case D: case N:ref0+=ce.getLength();break;
								case I: 
									{
									readseq.append(StringUtils.repeat(ce.getLength(),this.insert_base));
									// read0 += ce.getLength(); NO NO NO insert sequence is not provided in pslX format.
									break;
									}
								case M: case EQ: case X:
									{
									if(readBases==null)
										{
										readseq.append(genomicSeq.subSequence(ref0, ref0+ce.getLength())).toString().toUpperCase();
										}
									else //pslx
										{
										for(int j=0;j< ce.getLength();++j) {
											// I cannot understand why sometimes, there are fewer bases in the pslx format.
											if(read0+j < readBases.length())
												{
												readseq.append(readBases.charAt(read0+j));
												}
											else
												{
												readseq.append('N');
												}
											}
										}
									ref0+=ce.getLength();
									read0+=ce.getLength();
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
					rec.setAttribute("PG", prg.getId());
					rec.setCigar(new Cigar(cigar));
					sw.addAlignment(rec);
					}
			failed.close();failed=null;
			sw.close();sw=null;
			br.close();br=null;
			referenceSequenceFile.close();
			referenceSequenceFile=null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			err.printStackTrace();
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
