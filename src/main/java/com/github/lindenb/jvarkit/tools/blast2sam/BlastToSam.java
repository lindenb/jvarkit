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
package com.github.lindenb.jvarkit.tools.blast2sam;

import gov.nih.nlm.ncbi.blast.Hit;
import gov.nih.nlm.ncbi.blast.Hsp;
import gov.nih.nlm.ncbi.blast.Iteration;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.blast.BlastHspAlignment;



/**

BEGIN_DOC





### Example

The following Makefile downloads a reference , generates some FASTQs, align them with blastn and convert it to SAM:


```
BLASTN=/commun/data/packages/ncbi/ncbi-blast-2.2.28+/bin/blastn
SAMTOOLS=/commun/data/packages/samtools-0.1.19
JVARKIT=/home/lindenb/src/jvarkit-git/dist/
SHELL=/bin/bash
.PHONY:all reads clean
all: out.sam



out.sam: ref.fa ref.fa.fai out.read1.fq out.read2.fq
	paste \
		<(cat out.read1.fq | paste - - - - | cut -f 1,2 ) \
		<(cat out.read2.fq | paste - - - - | cut -f 1,2 ) |\
	tr "\t" "\n" |\
	sed 's/^@/>/' |\
	${BLASTN} -subject ref.fa -dust no -outfmt 5 | \
	java -jar ${JVARKIT}/blast2sam.jar -r ref.fa -p 500  |\
	${SAMTOOLS}/samtools view -Sh -f 2 - > $@
	
reads: out.read1.fq out.read2.fq
out.read1.fq out.read2.fq: ref.fa ref.fa.fai
	${SAMTOOLS}/misc/wgsim  -d 100 -N 500 -1 50 -2 50   $< out.read1.fq out.read2.fq > /dev/null

ref.fa:
	curl -k -o $@ "https://raw.github.com/lindenb/genomehub/master/data/rotavirus/rf/rf.fa"

ref.fa.fai: ref.fa
	${SAMTOOLS}/samtools faidx $<

clean:
	rm -f ref.fa.fai ref.fa out.sam 

```





### Output



```
@HD	VN:1.4	SO:unsorted
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:g1	LB:blast	DS:blast	SM:blast
@PG	ID:0	PN:blastn	VN:BLASTN_2.2.28+
@PG	ID:1	PN:com.github.lindenb.jvarkit.tools.blast2sam.BlastToSam	PP:0	VN:3365d9b714aa43d4fba44bfbf102a179a1f1573fCL:-r ref.fa -p 500
RF01_445_573_0:0:0_0:0:0_0/1	83	RF01	524	40	50=	=	445	-30	GTGCCTTGGTACACCATATTTATTTACTGTTGAAGCTACTATAGTGAATA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_445_573_0:0:0_0:0:0_0/2	163	RF01	445	40	50=	=	524	30	AATGCAGTTATGTTCTGGTTGGAAAAACATGAAAATGACGTTGCTGAAAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_1193_1294_1:0:0_1:0:0_1/1	83	RF01	1245	40	38=1X11=	=	1193	-3	CCATTACATGCATATTCTTTTTAGTCGAAAAAATTGTCATTCTACCAAAT	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_1193_1294_1:0:0_1:0:0_1/2	163	RF01	1193	40	4=1X45=	=	1245	3	CTGGATTACTATCAATGTCATCAGCGTCGAATGGTGAATCAAGACAACTA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_638_718_1:0:0_0:0:0_2/1	83	RF01	669	40	50=	=	638	18	ATGACAGTACTATCAGTTCTCTCGCAATTAAATAATCTTCATGAGAAAAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_638_718_1:0:0_0:0:0_2/2	163	RF01	638	40	4=1X45=	=	669	-18	CAAAATCTTCAATTGAAATGCTGATGTCAGTTTTTTCTCATGAAGATTAT	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_1404_1584_0:0:0_2:0:0_3/1	99	RF01	1404	40	50=	=	1535	179	ATTTATCTTACCATATGAATATTTCATAGCACAACATGCTGTAGTTGAAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_1404_1584_0:0:0_2:0:0_3/2	147	RF01	1535	40	1S42=1X6=	=	1404	-179	NGACACGTCTGTATATAGTACCATAGAGTTATTAGATAAAAAGGGTGTAA	#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:86.0662	BE:f:1.62562e-21	RG:Z:g1	NM:i:0	BS:f:46
RF01_284_373_0:0:0_1:0:0_5/1	99	RF01	284	40	50=	=	324	89	TAGTAAAATATGCAAAAGGTAAGCCGCTAGAAGCAGATTTGACAGTGAAT	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_284_373_0:0:0_1:0:0_5/2	147	RF01	324	40	8=1X41=	=	284	-89	AAAGTTCATATGTTATCTTGTTATTTTCATAATCCAACTCATTCACTGTC	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_1704_1823_1:0:0_0:0:0_7/1	83	RF01	1774	40	50=	=	1704	-21	ATTGAATTCGCTGCTTTCGTCTGCTTCTCTCCTGACGCTACAGCCCCATA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_1704_1823_1:0:0_0:0:0_7/2	163	RF01	1704	40	5=1X44=	=	1774	21	ACAGAGGCAAATTAATCTAATGGATTCATACGTTCAAATACCAGATGGTA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_689_741_1:0:0_1:0:0_8/1	83	RF01	692	40	19=1X30=	=	689	46	TGCCAGAGTCGATCTATTATAATATGACAGTACTATCAGTTCTCTCGCAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_689_741_1:0:0_1:0:0_8/2	163	RF01	689	40	30=1X19=	=	692	-46	TAATTGCGAGAGAACTGATAGTACTGTCATCTTCTAATAGATCGACTCTG	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_532_688_0:0:0_1:0:0_9/1	99	RF01	532	40	50=	=	639	156	ATAGTAGCTTCAACAGTAAATAAATATGGTGTACCAAGGCACAACGCGAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
(...)

```





END_DOC
*/


import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/*
BEGIN_DOC 

### Example

The following Makefile downloads a reference , generates some FASTQs, align them with blastn and convert it to SAM:

```
BLASTN=/commun/data/packages/ncbi/ncbi-blast-2.2.28+/bin/blastn
SAMTOOLS=/commun/data/packages/samtools-0.1.19
JVARKIT=/home/lindenb/src/jvarkit-git/dist/
SHELL=/bin/bash
.PHONY:all reads clean
all: out.sam



out.sam: ref.fa ref.fa.fai out.read1.fq out.read2.fq
	paste \
		<(cat out.read1.fq | paste - - - - | cut -f 1,2 ) \
		<(cat out.read2.fq | paste - - - - | cut -f 1,2 ) |\
	tr "\t" "\n" |\
	sed 's/^@/>/' |\
	${BLASTN} -subject ref.fa -dust no -outfmt 5 | \
	java -jar ${JVARKIT}/blast2sam.jar -r ref.fa -p 500  |\
	${SAMTOOLS}/samtools view -Sh -f 2 - > $@
	
reads: out.read1.fq out.read2.fq
out.read1.fq out.read2.fq: ref.fa ref.fa.fai
	${SAMTOOLS}/misc/wgsim  -d 100 -N 500 -1 50 -2 50   $< out.read1.fq out.read2.fq > /dev/null

ref.fa:
	curl -k -o $@ "https://raw.github.com/lindenb/genomehub/master/data/rotavirus/rf/rf.fa"

ref.fa.fai: ref.fa
	${SAMTOOLS}/samtools faidx $<

clean:
	rm -f ref.fa.fai ref.fa out.sam 

```



### Output

```
@HD	VN:1.4	SO:unsorted
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:g1	LB:blast	DS:blast	SM:blast
@PG	ID:0	PN:blastn	VN:BLASTN_2.2.28+
@PG	ID:1	PN:com.github.lindenb.jvarkit.tools.blast2sam.BlastToSam	PP:0	VN:3365d9b714aa43d4fba44bfbf102a179a1f1573f	CL:-r ref.fa -p 500
RF01_445_573_0:0:0_0:0:0_0/1	83	RF01	524	40	50=	=	445	-30	GTGCCTTGGTACACCATATTTATTTACTGTTGAAGCTACTATAGTGAATA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_445_573_0:0:0_0:0:0_0/2	163	RF01	445	40	50=	=	524	30	AATGCAGTTATGTTCTGGTTGGAAAAACATGAAAATGACGTTGCTGAAAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_1193_1294_1:0:0_1:0:0_1/1	83	RF01	1245	40	38=1X11=	=	1193	-3	CCATTACATGCATATTCTTTTTAGTCGAAAAAATTGTCATTCTACCAAAT	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_1193_1294_1:0:0_1:0:0_1/2	163	RF01	1193	40	4=1X45=	=	1245	3	CTGGATTACTATCAATGTCATCAGCGTCGAATGGTGAATCAAGACAACTA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_638_718_1:0:0_0:0:0_2/1	83	RF01	669	40	50=	=	638	18	ATGACAGTACTATCAGTTCTCTCGCAATTAAATAATCTTCATGAGAAAAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_638_718_1:0:0_0:0:0_2/2	163	RF01	638	40	4=1X45=	=	669	-18	CAAAATCTTCAATTGAAATGCTGATGTCAGTTTTTTCTCATGAAGATTAT	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_1404_1584_0:0:0_2:0:0_3/1	99	RF01	1404	40	50=	=	1535	179	ATTTATCTTACCATATGAATATTTCATAGCACAACATGCTGTAGTTGAAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_1404_1584_0:0:0_2:0:0_3/2	147	RF01	1535	40	1S42=1X6=	=	1404	-179	NGACACGTCTGTATATAGTACCATAGAGTTATTAGATAAAAAGGGTGTAA	#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:86.0662	BE:f:1.62562e-21	RG:Z:g1	NM:i:0	BS:f:46
RF01_284_373_0:0:0_1:0:0_5/1	99	RF01	284	40	50=	=	324	89	TAGTAAAATATGCAAAAGGTAAGCCGCTAGAAGCAGATTTGACAGTGAAT	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_284_373_0:0:0_1:0:0_5/2	147	RF01	324	40	8=1X41=	=	284	-89	AAAGTTCATATGTTATCTTGTTATTTTCATAATCCAACTCATTCACTGTC	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_1704_1823_1:0:0_0:0:0_7/1	83	RF01	1774	40	50=	=	1704	-21	ATTGAATTCGCTGCTTTCGTCTGCTTCTCTCCTGACGCTACAGCCCCATA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
RF01_1704_1823_1:0:0_0:0:0_7/2	163	RF01	1704	40	5=1X44=	=	1774	21	ACAGAGGCAAATTAATCTAATGGATTCATACGTTCAAATACCAGATGGTA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_689_741_1:0:0_1:0:0_8/1	83	RF01	692	40	19=1X30=	=	689	46	TGCCAGAGTCGATCTATTATAATATGACAGTACTATCAGTTCTCTCGCAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_689_741_1:0:0_1:0:0_8/2	163	RF01	689	40	30=1X19=	=	692	-46	TAATTGCGAGAGAACTGATAGTACTGTCATCTTCTAATAGATCGACTCTG	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:87.9128	BE:f:4.51982e-22	RG:Z:g1	NM:i:0	BS:f:47
RF01_532_688_0:0:0_1:0:0_9/1	99	RF01	532	40	50=	=	639	156	ATAGTAGCTTCAACAGTAAATAAATATGGTGTACCAAGGCACAACGCGAA	JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ	BB:f:93.4528	BE:f:9.71473e-24	RG:Z:g1	NM:i:0	BS:f:50
(...)

```



END_DOC
*/
@Program(name="blast2sam",
description="Convert a **BLASTN-XML** input to SAM",
keywords={"sam","blast"})
public class BlastToSam extends Launcher
	{
	private static final Logger LOG = Logger.build(BlastToSam.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-p","--expect_size"},description="input is an interleaved list of sequences forward and reverse (paired-ends). 0: not interleaved")
	private int EXPECTED_SIZE = 0 ;

	@Parameter(names={"-r","--reference"},description="Indexed fasta Reference")
	private File faidx = null;
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	

	private SAMSequenceDictionary dictionary;
	private Unmarshaller unmarshaller;
	//fool javac
	@SuppressWarnings("unused")
	private final static gov.nih.nlm.ncbi.blast.ObjectFactory _foolJavac=null;
	
	private static class SequenceIteration
		{
		//Iteration iteration;
		//String queryDef;
		List<SAMRecord> records=new ArrayList<SAMRecord>();
		}
	
	private class Paired
		implements Comparable<Paired>
		{
		SAMRecord rec1;
		SAMRecord rec2;
		
		void completeFlags()
			{
			rec1.setFirstOfPairFlag(true);
			rec1.setSecondOfPairFlag(false);
			
			rec2.setFirstOfPairFlag(false);
			rec2.setSecondOfPairFlag(true);
			

			rec1.setReadPairedFlag(true);
			rec2.setReadPairedFlag(true);
			
			rec1.setProperPairFlag(false);
			rec2.setProperPairFlag(false);
			
			for(int i=0;i< 2;++i)
				{
				SAMRecord recA=(i==0?rec1:rec2);
				if(recA.getReadUnmappedFlag())
					{
					recA.setReferenceIndex(-1);
					}
				}
			
			for(int i=0;i< 2;++i)
				{
				SAMRecord recA=(i==0?rec1:rec2);
				SAMRecord recB=(i==0?rec2:rec1);
				recA.setMateUnmappedFlag(recB.getReadUnmappedFlag());
				
				
				if(!recB.getReadUnmappedFlag())
					{
					recA.setMateReferenceName(recB.getReferenceName());
					recA.setMateReferenceIndex(recB.getReferenceIndex());
					recA.setMateNegativeStrandFlag(recB.getReadNegativeStrandFlag());
					recA.setMateAlignmentStart(recB.getAlignmentStart());
					}
				else
					{
					recA.setMateReferenceIndex(-1);
					
					}
				}	
			if( !rec1.getReadUnmappedFlag() &&
				!rec2.getReadUnmappedFlag() && 
				rec1.getReadNegativeStrandFlag()!=rec2.getReadNegativeStrandFlag() &&
				rec1.getReferenceName().equals(rec2.getReferenceName())		 
				)
				{
				int len=rec1.getAlignmentStart()-rec2.getAlignmentEnd();
				if(Math.abs(len) <= EXPECTED_SIZE)
					{
					rec1.setProperPairFlag(true);
					rec2.setProperPairFlag(true);
					}
				rec1.setInferredInsertSize(-len);
				rec2.setInferredInsertSize(len);
				
				}
			}
		
		int score()
			{
			int v=0;
			v+=rec1.getReadUnmappedFlag()?0:1;
			v+=rec2.getReadUnmappedFlag()?0:1;
			if( !rec1.getReadUnmappedFlag() && !rec2.getReadUnmappedFlag()  )		 
				{
				if(rec1.getReferenceName().equals(rec2.getReferenceName()))
					{
					v+=2;
					if(rec1.getReadNegativeStrandFlag()!=rec2.getReadNegativeStrandFlag())
						{
						v+=4;
						int len=rec1.getAlignmentStart()-rec2.getAlignmentEnd();
						if(Math.abs(len) <= EXPECTED_SIZE)
							{
							v+=8;
							}
						
						}
					}
				}
			return v;
			}
		
		@Override
		public int compareTo(Paired p)
			{
			return p.score()-score();//greater is best
			}
		}
	
	private BlastToSam()
		{
		
		}
	
	
	

	private Iteration peekIteration(XMLEventReader r) throws XMLStreamException,JAXBException
		{
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
		
			if(!(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals("Iteration")))
				{
				r.next();
				continue;
				}
			return this.unmarshaller.unmarshal(r, Iteration.class).getValue();
			}
		return null;
		}
	
	private void fillHeader(XMLEventReader r,SAMProgramRecord prog) throws XMLStreamException,JAXBException
		{
		while(r.hasNext())
			{
			XMLEvent evt=r.peek();
		
			if(!(evt.isStartElement()))
				{
				r.next();
				continue;
				}
			StartElement E=evt.asStartElement();
			String name=E.getName().getLocalPart();
			if(name.equals("BlastOutput_iterations")) break;
			r.next();
			if(name.equals("BlastOutput_program"))
				{
				prog.setProgramName(r.getElementText());
				}
			else if(name.equals("BlastOutput_version"))
				{
				prog.setProgramVersion(r.getElementText().replace(' ', '_'));
				}
			}
		}
	
	private void dumpSingle(final SAMFileWriter w,final SequenceIteration si)
		{
		boolean first=true;
		for(final SAMRecord rec:si.records)
			{
			rec.setSecondaryAlignment(!first);
			first=false;
			w.addAlignment(rec);
			}
		si.records.clear();
		}
	
	private void run_single(
			final SAMFileWriter w,
			final XMLEventReader r,
			final SAMFileHeader header
			)
			throws XMLStreamException,JAXBException
		{
		final List<Iteration> stack=new ArrayList<Iteration>();
		String prev=null;
		for(;;)
			{
			final Iteration iter1=peekIteration(r);
			if(iter1==null || !(iter1.getIterationQueryDef().equals(prev)))
				{
				final SequenceIteration si=convertIterationToSequenceIteration(stack,header);
				dumpSingle(w,si);
				if(iter1==null) break;
				stack.clear();
				prev=iter1.getIterationQueryDef();
				}
			stack.add(iter1);
			}
		}
	
	private SequenceIteration convertIterationToSequenceIteration(
			final List<Iteration> stack,
			final SAMFileHeader header
			)
			throws XMLStreamException,JAXBException
			{
			final SequenceIteration sequenceIteration=new SequenceIteration(); 
			if(stack.isEmpty()) return sequenceIteration;
			
			final SAMReadGroupRecord rg1=header.getReadGroup("g1");
			//sequenceIteration.iteration=iter1;
			
			final SAMRecordFactory samRecordFactory=new DefaultSAMRecordFactory();

			
			
			final StringBuilder readContent=new StringBuilder();
			final int iterLength=Integer.parseInt(stack.get(0).getIterationQueryLen());
			
			for(final Iteration iter1:stack)
				{
				for(final Hit hit: iter1.getIterationHits().getHit())
					{
					for(final Hsp hsp: hit.getHitHsps().getHsp())
						{
						for(final BlastHspAlignment.Align a:new BlastHspAlignment(hsp))
							{
							char c=a.getQueryChar();
							if(!Character.isLetter(c)) continue;
							final int queryIndex0=a.getQueryIndex1()-1;
							while(readContent.length()<=queryIndex0) readContent.append('N');
							if(readContent.charAt(queryIndex0)=='N')
								{
								readContent.setCharAt(queryIndex0, c);
								}
							else if(readContent.charAt(queryIndex0)!=c)
								{
								throw new IllegalStateException(
									"Expected character '"+readContent.charAt(queryIndex0)+"' but got '"+c+"' at "+queryIndex0+"\n"+
									hsp.getHspQseq()+"\n"+
									hsp.getHspMidline()+"\n"+
									hsp.getHspHseq()+"\n"+
									readContent+"\n"
									);
								}
							}
						}
					}
				}
			
			
			for(Iteration iter1:stack)
				{
				for(Hit hit: iter1.getIterationHits().getHit())
					{
					for(Hsp hsp: hit.getHitHsps().getHsp())
						{
						SAMRecord rec=samRecordFactory.createSAMRecord(header);
						rec.setReadUnmappedFlag(false);
						rec.setReadName(iter1.getIterationQueryDef());
						if( hit.getHitAccession()!=null &&
							!hit.getHitAccession().trim().isEmpty() &&
							this.dictionary.getSequence(hit.getHitAccession())!=null
							)
							{
							rec.setReferenceName(hit.getHitAccession());
							}
						else
							{
							rec.setReferenceName(hit.getHitDef());
							}
						final SAMSequenceRecord ssr=this.dictionary.getSequence(hit.getHitDef());
						if(ssr==null)
							{
							LOG.warn("Hit is not in SAMDictionary "+hit.getHitDef());
							rec.setReferenceIndex(-1);
							}
						else
							{
							rec.setReferenceIndex(ssr.getSequenceIndex());
							}
						
						final BlastHspAlignment blastHspAlignment=new BlastHspAlignment(hsp);
						rec.setReadNegativeStrandFlag(blastHspAlignment.isPlusMinus());
	
						
						final List<CigarOperator> cigarL=new ArrayList<CigarOperator>();
						for(BlastHspAlignment.Align a:blastHspAlignment)
							{
							//System.err.println("##"+a);
							if(a.getMidChar()=='|')
								{
								cigarL.add(CigarOperator.EQ);
								}
							else if(a.getMidChar()==':')
								{
								cigarL.add(CigarOperator.M);
								}
							else if(a.getHitChar()=='-')
								{
								cigarL.add(CigarOperator.I);
								}
							else if(a.getQueryChar()=='-')
								{
								cigarL.add(CigarOperator.D);
								}
							else
								{
								cigarL.add(CigarOperator.X);
								}
	
							}
	
						
						if(cigarL.size()!=hsp.getHspMidline().length())
							{
							throw new IllegalStateException("Boumm");
							}
						
						
						Cigar cigarE=new Cigar();
						
						if(blastHspAlignment.getQueryFrom1()>1)
							{
							cigarE.add(new CigarElement(
									blastHspAlignment.getQueryFrom1()-1,
									CigarOperator.S
									));
							}
						int x=0;
						while(x< cigarL.size())
							{
							int y=x+1;
							while(y< cigarL.size() && cigarL.get(x)==cigarL.get(y))
								{
								++y;
								}
							cigarE.add(new CigarElement(y-x, cigarL.get(x)));
							x=y;
							}
						/* soft clip */ 
						if(blastHspAlignment.getQueryTo1()< readContent.length())
							{
							cigarE.add(new CigarElement(
									(readContent.length()-blastHspAlignment.getQueryTo1()),
									CigarOperator.S 
									));
							}
						/* hard clip */
						if(readContent.length() < iterLength)
							{
							cigarE.add(new CigarElement(
									(iterLength-readContent.length()),
									CigarOperator.H
									));
							}
						
						
						rec.setCigar(cigarE);
						rec.setMappingQuality(40);
						rec.setAlignmentStart(Math.min(blastHspAlignment.getHitFrom1(),blastHspAlignment.getHitTo1()));
						rec.setAttribute("BB", Float.parseFloat(hsp.getHspBitScore()));
						rec.setAttribute("BE", Float.parseFloat(hsp.getHspEvalue()));
						rec.setAttribute("BS", Float.parseFloat(hsp.getHspScore()));
						rec.setAttribute("NM", Integer.parseInt(hsp.getHspGaps()));
						rec.setAttribute("RG", rg1.getId());
						// setAlignmentEnd not supported in SAM API
						//rec.setAlignmentEnd(Math.max(blastHspAlignment.getHitFrom1(),blastHspAlignment.getHitTo1())); 
						sequenceIteration.records.add(rec);
						}
					}
				}
			
			if(readContent.length()==0)
				{
				readContent.append('N');
				}
			
			byte readBases[]=readContent.toString().getBytes();
			char readQuals[]=new char[readBases.length];
			
			 
			
			for(int i=0;i< readBases.length;++i)
				{
				readQuals[i]=(readBases[i]=='N'?'#':'J');
				}

			
			
			if(sequenceIteration.records.isEmpty())
				{
				SAMRecord rec=samRecordFactory.createSAMRecord(header);
				rec.setReadName(stack.get(0).getIterationQueryDef());
				rec.setReadUnmappedFlag(true);
				rec.setAttribute("RG", rg1.getId());
				sequenceIteration.records.add(rec);
				}
			
			
				
			
			for(SAMRecord rec:sequenceIteration.records)
				{
				rec.setReadString(new String(readBases));
				rec.setReadBases(readBases);
				rec.setBaseQualityString(new String(readQuals,0,readQuals.length));
				rec.setBaseQualities(htsjdk.samtools.SAMUtils.fastqToPhred(new String(readQuals,0,readQuals.length)));
				}
		return sequenceIteration;
		}
	
	private static SAMRecord cloneSAMRecord(final SAMRecord rec)
		{
		try {
			return (SAMRecord)rec.clone();
			}
		catch (Exception e)
			{
			throw new RuntimeException("Cannot clone a SAMRecord ?",e);
			}
		}
	
	private void dumpPaired(SAMFileWriter w,SequenceIteration si1,SequenceIteration si2)
		{
		if(si1.records.isEmpty()) return;
		
		SequenceIteration siL[]=new SequenceIteration[]{si1,si2};
		for(SequenceIteration si:siL)
			{
			for(SAMRecord rec:si.records)
				{
				rec.setReadPairedFlag(true);
				rec.setMateUnmappedFlag(true);
				}
			}
		List<Paired> paired=new ArrayList<Paired>();
		int x=0;
		while(x < si1.records.size())
			{
			SAMRecord rec1=si1.records.get(x);
			int y=0;
			
			while(y < si2.records.size())
				{
				SAMRecord rec2=si2.records.get(y);
				
				Paired pair=new Paired();
				
				pair.rec1=cloneSAMRecord(rec1);
				pair.rec2=cloneSAMRecord(rec2);
				
				if(pair.rec1.getReadUnmappedFlag() && pair.rec2.getReadUnmappedFlag())
					{
					++y;
					continue;
					}	
				
				if(!paired.isEmpty() )
					{
					int cmp= pair.compareTo(paired.get(0)) ;
					if(cmp<0)
						{
						paired.clear();
						}
					else if(cmp>0)
						{
						++y;
						continue;
						}
					}
				paired.add(pair);
				++y;
				}
			++x;
			}
		
		if(paired.isEmpty())
			{
			Paired pair=new Paired();
			pair.rec1=cloneSAMRecord(si1.records.get(0));
			pair.rec2=cloneSAMRecord(si2.records.get(0));
			paired.add(pair);
			}
		
		for(int i=0;i< paired.size();++i)
			{
			Paired pair=paired.get(i);
			pair.completeFlags();
			if(!pair.rec1.getReadUnmappedFlag())
				{
				pair.rec1.setSecondaryAlignment(i!=0);
				}
			if(!pair.rec2.getReadUnmappedFlag())
				{
				pair.rec2.setSecondaryAlignment(i!=0);
				}
			w.addAlignment(pair.rec1);
			w.addAlignment(pair.rec2);
			}
		si1.records.clear();
		si2.records.clear();
		}
	
	private void run_paired(
			SAMFileWriter w,
			XMLEventReader r,
			SAMFileHeader header
			)
			throws XMLStreamException,JAXBException
		{
		List<Iteration> stack1=new ArrayList<Iteration>();
		Iteration iter=null;
		for(;;)
			{
			String prev_name=null;
			if( iter==null)
				{
				iter=peekIteration(r);
				if(iter==null) break;
				}
			stack1.add(iter);
			List<Iteration> stack2=new ArrayList<Iteration>();
			prev_name=iter.getIterationQueryDef();
			
			//pileup first of pair
			for(;;)
				{
				iter=peekIteration(r);
				if(iter==null)
					{
					throw new RuntimeException("Illegal number of read forward/reverse");
					}
				else if(iter.getIterationQueryDef().equals(prev_name))
					{
					stack1.add(iter);
					}
				else
					{
					stack2.add(iter);
					prev_name=iter.getIterationQueryDef();
					break;
					}
				}
			
			//pileup second of pair
			for(;;)
				{
				iter=peekIteration(r);
				if(iter==null || !iter.getIterationQueryDef().equals(prev_name))
					{
					SequenceIteration si1=convertIterationToSequenceIteration(stack1, header);
					SequenceIteration si2=convertIterationToSequenceIteration(stack2, header);
					dumpPaired(w,si1,si2);
					stack1.clear();
					stack2.clear();
					break;
					}
				else
					{
					stack2.add(iter);
					}
				}
			if(iter==null) break;
			}
		
		}
	
	@Override
	public int doWork(List<String> args) {
		
		if(this.faidx==null || !this.faidx.exists() || !this.faidx.isFile()) {
			LOG.error("Option -R was not defined or dictionary missing");
			return -1;
		}
		final boolean interleaved_input=this.EXPECTED_SIZE>0;
		final int maxRecordsInRam=5000;
		SAMFileWriter sfw=null;
		XMLEventReader rx=null;
		final SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
		sfwf.setCreateIndex(false);
		sfwf.setMaxRecordsInRam(maxRecordsInRam);
		sfwf.setCreateMd5File(false);
		sfwf.setUseAsyncIo(false);
		final SAMFileHeader header=new SAMFileHeader();
		try
			{
			LOG.info("opening "+faidx);
			this.dictionary=SAMSequenceDictionaryExtractor.extractDictionary(faidx.toPath());
			header.setSortOrder(SortOrder.unsorted);
			header.setSequenceDictionary(this.dictionary);
			
			
			final JAXBContext jc = JAXBContext.newInstance("gov.nih.nlm.ncbi.blast");
			this.unmarshaller=jc.createUnmarshaller();
			final XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.FALSE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			xmlInputFactory.setXMLResolver(new XMLResolver()
				{
				@Override
				public Object resolveEntity(String arg0, String arg1, String arg2,
						String arg3) throws XMLStreamException
					{
					LOG.info("resolveEntity:" +arg0+"/"+arg1+"/"+arg2);
					return null;
					}
				});
			final String inputName=oneFileOrNull(args);
			if(inputName==null)
				{
				LOG.info("Reading from stdin");
				rx=xmlInputFactory.createXMLEventReader(stdin());
				}
			else if(args.size()==1)
				{
				LOG.info("Reading from "+inputName);
				rx=xmlInputFactory.createXMLEventReader(IOUtils.openURIForBufferedReading(inputName));
				}
			else
				{
				LOG.error("Illegal number of args");
				return -1;
				}
			
			
			final SAMProgramRecord prg2=header.createProgramRecord();
			fillHeader(rx,prg2);
			final SAMProgramRecord prg1=header.createProgramRecord();
			prg1.setCommandLine(getProgramCommandLine());
			prg1.setProgramVersion(getVersion());
			prg1.setProgramName(getProgramName());
			prg1.setPreviousProgramGroupId(prg2.getId());
			final SAMReadGroupRecord rg1=new SAMReadGroupRecord("g1");
			rg1.setLibrary("blast");
			rg1.setSample("blast");
			rg1.setDescription("blast");
			header.addReadGroup(rg1);
			
			sfw = this.writingBamArgs.openSAMFileWriter(outputFile,header, true);
			
			if(interleaved_input)
				{
				run_paired(sfw,rx,header);
				}
			else
				{
				run_single(sfw,rx,header);
				}
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}	
		finally
			{
			CloserUtil.close(sfw);
			CloserUtil.close(rx);
			}
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BlastToSam().instanceMainWithExit(args);
		}

	}
