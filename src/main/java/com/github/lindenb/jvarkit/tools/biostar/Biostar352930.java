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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ReadNameSortMethod;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

/** 
BEGIN_DOC


## Motivation

Fills the empty SEQ(*) and QUAL(*) in a bam file using the the reads with the same name carrying this information.

## Input

a sam or a bam file. Input **MUST** be sorted on query name using picard. ( see see https://github.com/samtools/hts-specs/issues/5 )

## Example

A read `A00223:8:H7YG3DMXX:1:1101:1163:35383` in the remote bam file below is missing SEQ/QUAL while another has this information

```
$ wget -q -O - "https://gist.githubusercontent.com/toddknutson/90430a0dd736898037ed18bcd044df7f/raw/87c1ea5a548ac71c628d2b72f5bd6ee6415efbcd/gistfile1.txt" | grep -F "A00223:8:H7YG3DMXX:1:1101:1163:35383"
A00223:8:H7YG3DMXX:1:1101:1163:35383	16	chr7	6027025	9	151M	*	0	0	GCTAGAAGACAGCAGACCCCTTGTCTGTCCTAGAGGGCTCCTTCTTGGTTCTGGAGTCTTTGGGCTGTGAGGCTTGTTCTCTGTTGTGTGACGAAGAGAAAAGGCCTCTCGCAGTCTGGAAATGGACACGTCTTTTTTTTCTTCTCCAGTC	FFFFFFFFFFFFFF,:FFFF,F,FFFFFFFF:FFFFFFFFF:FF::FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFF	NM:i:2	MD:Z:14T7T128	AS:i:2978	XS:i:2894	RG:Z:INV_4SMALL_Aug_17_2017_1
A00223:8:H7YG3DMXX:1:1101:1163:35383	256	chr7	6776783	0	151M	*	0	0	*	*	NM:i:6	MD:Z:17G0G109A7A2T0C10	AS:i:2894	RG:Z:INV_4SMALL_Aug_17_2017_1

```

get the sam/bam file, sort it on queryname using picard

```
$ wget -q -O - "https://gist.githubusercontent.com/toddknutson/90430a0dd736898037ed18bcd044df7f/raw/87c1ea5a548ac71c628d2b72f5bd6ee6415efbcd/gistfile1.txt" |\
 	java -jar /path/to/picard.jar SortSam I=/dev/stdin O=/dev/stdout SO=queryname VALIDATION_STRINGENCY=LENIENT |\
 	java -jar dist/biostar352930.jar
```

output:

```
@HD	VN:1.5	SO:queryname
(...)
@RG	ID:INV_4SMALL_Aug_17_2017_1	SM:0X	PL:ILLUMINA	LB:libeX
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -t 12 -R @RG\tID:sample1\tSM:0X\tPL:ILLUMINA\tLB:libeX -Y -a -k 11 -O 2 -B 1 -A 20 hg19_canonical_PARmaskedonchrY.fasta R1.fastq
@CO	biostar352930. compilation: 2018-12-05:10-12-50 githash: 03f79caf6a62458937b2d55d0661737617334733 htsjdk: 2.15.0. cmd:
A00223:8:H7YG3DMXX:1:1101:1163:35383	256	chr7	6776783	0	151M	*	0	0	GACTGGAGAAGAAAAAAAAGACGTGTCCATTTCCAGACTGCGAGAGGCCTTTTCTCTTCGTCACACAACAGAGAACAAGCCTCACAGCCCAAAGACTCCAGAACCAAGAAGGAGCCCTCTAGGACAGACAAGGGGTCTGCTGTCTTCTAGC	FFFFFF,FF:FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF::FF:FFFFFFFFF:FFFFFFFF,F,FFFF:,FFFFFFFFFFFFFF	MD:Z:17G0G109A7A2T0C10	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:6	AS:i:2894
A00223:8:H7YG3DMXX:1:1101:1163:35383	16	chr7	6027025	9	151M	*	0	0	GCTAGAAGACAGCAGACCCCTTGTCTGTCCTAGAGGGCTCCTTCTTGGTTCTGGAGTCTTTGGGCTGTGAGGCTTGTTCTCTGTTGTGTGACGAAGAGAAAAGGCCTCTCGCAGTCTGGAAATGGACACGTCTTTTTTTTCTTCTCCAGTC	FFFFFFFFFFFFFF,:FFFF,F,FFFFFFFF:FFFFFFFFF:FF::FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF:FF,FFFFFF	MD:Z:14T7T128	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:2	AS:i:2978	XS:i:2894
A00223:8:H7YG3DMXX:1:1101:15338:10113	256	chr7	6777077	0	151M	*	0	0	GTTCAGCATCCCAGACACGGGCAGTCACTGCAGCAGCGAGTATGCGGCCAGCTCCCCAGGGGACAGGGGCTCGCAGGAACATGTGGACTCTCAGGAGAAAGCGCCTGAAACTGACGACTCTTTTTCAGATGTGGACTGCCATTCAAACCAG	,FFFFF,F,FFFFFFFFFFFF:F,F:F::F,FFFFFFFFFFFFFFFFFFFFFFF,F:FFF,:FFFFFFFFFF,FFF:FFFFF,FF,F:FFFFFF:F,FFFF:,::,,FFFFF:FFFF::F,FFF,:,,F:,FF,F,FF,FFF,FFFF,FF,	MD:Z:41G2T7A98	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:3	AS:i:2957
A00223:8:H7YG3DMXX:1:1101:15338:10113	16	chr7	6026731	4	151M	*	0	0	CTGGTTTGAATGGCAGTCCACATCTGAAAAAGAGTCGTCAGTTTCAGGCGCTTTCTCCTGAGAGTCCACATGTTCCTGCGAGCCCCTGTCCCCTGGGGAGCTGGCCGCATACTCGCTGCTGCAGTGACTGCCCGTGTCTGGGATGCTGAAC	,FF,FFFF,FFF,FF,F,FF,:F,,:,FFF,F::FFFF:FFFFF,,::,:FFFF,F:FFFFFF:F,FF,FFFFF:FFF,FFFFFFFFFF:,FFF:F,FFFFFFFFFFFFFFFFFFFFFFF,F::F:F,F:FFFFFFFFFFFF,F,FFFFF,	MD:Z:44T106	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:1	AS:i:2999	XS:i:2957
A00223:8:H7YG3DMXX:1:1101:15881:12242	0	chr7	6026908	11	151M	*	0	0	GTGCCCCGAGTCCTTCTCCACCTCCGCTCTGTCCGTAGGGTCACTGGGTCCGTGACTGGAACTCACTGCCTCTTTCTGAGATCTCAGGACGCCTTTGTCAGAGATGGCACCTGAAGTGCTAGAAGACAGCATACCCCTTTTCTGTCCTAGA	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,F:FFFFFFFFFFF:F,FFFFFFFFFFFFFFF:FFFFFFFFF:FF:FFF:FFFF:FFFFFFFF,FFFFFFFFF	MD:Z:80G70	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:1	AS:i:2999	XS:i:2895
A00223:8:H7YG3DMXX:1:1101:15881:12242	272	chr7	6776900	0	151M	*	0	0	TCTAGGACAGAAAAGGGGTATGCTGTCTTCTAGCACTTCAGGTGCCATCTCTGACAAAGGCGTCCTGAGATCTCAGAAAGAGGCAGTGAGTTCCAGTCACGGACCCAGTGACCCTACGGACAGAGCGGAGGTGGAGAAGGACTCGGGGCAC	FFFFFFFFF,FFFFFFFF:FFFF:FFF:FF:FFFFFFFFF:FFFFFFFFFFFFFFF,F:FFFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	MD:Z:22T0C17A28C28G50T0	RG:Z:INV_4SMALL_Aug_17_2017_1	NM:i:6	AS:i:2895
```

the read `A00223:8:H7YG3DMXX:1:1101:1163:35383` is now fixed.


END_DOC

*/

@Program(name="biostar352930",
	description="Fills the empty SEQ(*) and QUAL(*) in a bam file using the the reads with the same name carrying this information.",
	keywords= {"sam","bam"},
	biostars=352930
	)
public class Biostar352930 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar352930.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	
	@Override
	public int doWork(final List<String> args) {
		SamReader samReader = null;
		SAMFileWriter samFileWriter = null;
		try {
			long count_rescued_seq = 0L;
			long count_rescued_qual = 0L;
			final SamReaderFactory srf = super.createSamReaderFactory().
						validationStringency(ValidationStringency.SILENT);
			final String input = oneFileOrNull(args);
			if(input==null)
				{
				samReader =srf.open(SamInputResource.of(stdin()));
				}
			else
				{
				samReader =srf.open(SamInputResource.of(input));
				}
			final SAMFileHeader header = samReader.getFileHeader();
			if(header.getSortOrder()!=SAMFileHeader.SortOrder.queryname) {
				LOG.error("Expected SAM input to be sorted on "+SAMFileHeader.SortOrder.queryname+" but got "+header.getSortOrder() );
				return -1;
				}
			JVarkitVersion.getInstance().addMetaData(this, header);
			samFileWriter = this.writingBamArgs.openSAMFileWriter(this.outputFile, header, true);
			final CloseableIterator<SAMRecord> iter = samReader.iterator();
			final Comparator<SAMRecord> readNameCmp = ReadNameSortMethod.picard.get();
			final Iterator<List<SAMRecord>> eq_range =  new EqualRangeIterator<>(iter,readNameCmp);
			
			final Predicate<SAMRecord> noCigarOrNoHardClip = R->R.getCigar()==null || R.getCigar().getCigarElements().stream().noneMatch(C->C.getOperator().equals(CigarOperator.H));
			
			while(eq_range.hasNext()) {
				final List<SAMRecord> records = eq_range.next();
				
				for(int choice=0;choice<3;++choice) {
					final Predicate<SAMRecord> selectRead;
					switch(choice)
						{	
						case 0: selectRead = R->!R.getReadPairedFlag();break;
						case 1: selectRead = R->R.getReadPairedFlag() && R.getFirstOfPairFlag();break;
						case 2: selectRead = R->R.getReadPairedFlag() && R.getSecondOfPairFlag();break;
						default: throw new IllegalStateException();
						}
					
					// get the strand+ SEQUENCE string
					final String unclippedseq = records.
						stream().
						filter(selectRead).
						filter(noCigarOrNoHardClip).
						filter(R->!R.getReadString().equals(SAMRecord.NULL_SEQUENCE_STRING)).
						map(R->{
							if(R.getReadNegativeStrandFlag())
								{
								return  SequenceUtil.reverseComplement(R.getReadString());
								}
							else
								{
								return R.getReadString();
								}
							}).
						findFirst().
						orElse(null);
					
					// get the strand+ QUAL string
					final String unclippedqual = records.
						stream().
						filter(selectRead).
						filter(noCigarOrNoHardClip).
						filter(R->!R.getBaseQualityString().equals(SAMRecord.NULL_QUALS_STRING)).
						map(R->{
							if(R.getReadNegativeStrandFlag())
								{
								return  new StringBuilder(R.getBaseQualityString()).reverse().toString();
								}
							else
								{
								return R.getBaseQualityString();
								}
							}).
						findFirst().
						orElse(null);
					
					// loop over the record and fix
					for(final SAMRecord rec: records)
						{
						if(!selectRead.test(rec)) continue;
						
						String clippedseq = unclippedseq;
						String clippedqual = unclippedqual;
						
						if(rec.getReadUnmappedFlag() && rec.getCigar()!=null &&
							rec.getCigar().numCigarElements()>1)
							{
							final Cigar cigar = rec.getCigar();
							CigarElement ce = cigar.getLastCigarElement();
							if(ce.getOperator().equals(CigarOperator.H)) {
								clippedseq = clippedseq.substring(0,clippedseq.length()-ce.getLength());
								clippedqual = clippedqual.substring(0,clippedqual.length()-ce.getLength());
								}
							ce = cigar.getFirstCigarElement();
							if(ce.getOperator().equals(CigarOperator.H)) {
								clippedseq = clippedseq.substring(ce.getLength());
								clippedqual = clippedqual.substring(ce.getLength());
								}
							}
						
						if(!StringUtil.isBlank(clippedseq) && 
							rec.getReadString().equals(SAMRecord.NULL_SEQUENCE_STRING))
							{
							if(rec.getReadNegativeStrandFlag())
								{
								rec.setReadString(SequenceUtil.reverseComplement(clippedseq));
								}
							else
								{
								rec.setReadString(clippedseq);
								}
							++count_rescued_seq;
							}
						if(!StringUtil.isBlank(clippedqual) && 
							rec.getBaseQualityString().equals(SAMRecord.NULL_QUALS_STRING))
							{
							if(rec.getReadNegativeStrandFlag())
								{
								rec.setBaseQualityString(new StringBuilder(clippedqual).reverse().toString());
								}
							else
								{
								rec.setBaseQualityString(clippedqual);
								}
							++count_rescued_qual;
							}
						}
					}
				for(final SAMRecord R: records) samFileWriter.addAlignment(R);
				}
			CloserUtil.close(eq_range);
 			iter.close();
			samFileWriter.close();samFileWriter=null;
			samReader.close();samReader=null;
			LOG.info("Done : rescued seq: "+count_rescued_seq +" rescued qual: "+count_rescued_qual );
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
		} finally
			{
			CloserUtil.close(samReader);
			CloserUtil.close(samFileWriter);
			}
	}
		
	
	public static void main(final String[] args) throws IOException
		{
		new Biostar352930().instanceMainWithExit(args);
		}
		

	}
