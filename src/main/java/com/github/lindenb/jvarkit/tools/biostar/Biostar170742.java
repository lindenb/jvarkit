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
package com.github.lindenb.jvarkit.tools.biostar;
import java.io.PrintStream;
import java.nio.file.Path;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.MultiBamLauncher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
BEGIN_DOC

## Example


```
$ java -jar dist-2.0.1/biostar170742.jar \
	-R ref.fa \
	S1.bam  | head -n 15
 
0 rotavirus 1 70 rotavirus_1_317_5:0:0_7:0:0_2de/1 + 60
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT

1 rotavirus 1 70 rotavirus_1_535_4:0:0_4:0:0_1a6/2 + 60
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT
GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT

2 rotavirus 1 70 rotavirus_1_543_5:0:0_11:0:0_390/2 + 60
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT
GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT

3 rotavirus 1 70 rotavirus_1_578_3:0:0_7:0:0_7c/1 + 60
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTATTATT
GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATT
(...)
```

END_DOC
 */
@Program(name="biostar170742",
description="convert sam format to axt Format",
biostars=170742,
modificationDate="20210412",
creationDate="20151228",
keywords={"sam","axt"})
public class Biostar170742 extends MultiBamLauncher {

	private static final Logger LOG = Logger.build(Biostar170742.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int processInput(final SAMFileHeader header, final CloseableIterator<SAMRecord> iter) {
		PrintStream out=null;
		GenomicSequence genomicSequence = null;
		ReferenceSequenceFile indexedFastaSequenceFile= null;
		try
			{
			indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(getRequiredReferencePath());
			long align_id=0;
			final StringBuilder refseq = new StringBuilder();
			final StringBuilder readseq = new StringBuilder();
			out =  super.openPathOrStdoutAsPrintStream(this.outputFile);
			
			while(iter.hasNext())
				{
				final SAMRecord rec= iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar==null) continue;
				final byte readbases[] = rec.getReadBases();
				if(readbases==null || readbases==SAMRecord.NULL_SEQUENCE) continue;

				if(genomicSequence==null || !rec.getReferenceName().equals(genomicSequence.getChrom()))
					{
					genomicSequence = new  GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
					}
				
				int refpos1 = rec.getAlignmentStart();
				int readpos = 0; 
				refseq.setLength(0);
				readseq.setLength(0);
				for(final CigarElement ce: cigar.getCigarElements())
					{
					final CigarOperator op = ce.getOperator();
					if(op.equals(CigarOperator.S))
						{
						readpos+=ce.getLength();
						continue;
						}
					if(op.equals(CigarOperator.H))
						{
						continue;
						}
				

					
					for(int i=0;i< ce.getLength();++i)
						{
						if( op.consumesReferenceBases() &&
							op.consumesReadBases())
							{
							refseq.append(genomicSequence.charAt(refpos1 - 1));
							readseq.append((char)readbases[readpos]);
							readpos++;
							refpos1++;
							}
						else if( op.consumesReferenceBases() )
							{
							refseq.append(genomicSequence.charAt(refpos1 -1));
							readseq.append('-');
							refpos1++;
							}
						else if( op.consumesReadBases() )
							{
							refseq.append('-');
							readseq.append((char)readbases[readpos]);
							readpos++;
							}
						}
					}
				out.print(align_id);
				out.print(' ');
				out.print(rec.getReferenceName());
				out.print(' ');
				out.print(rec.getAlignmentStart());
				out.print(' ');
				out.print(rec.getAlignmentEnd());
				out.print(' ');
				out.print(rec.getReadName());
				if(rec.getReadPairedFlag())
					{
					if(rec.getFirstOfPairFlag())
						{
						out.print("/1");	
						}
					else if(rec.getSecondOfPairFlag())
						{
						out.print("/2");	
						}
					}
				out.print(' ');
				out.print(1+rec.getAlignmentStart()-rec.getUnclippedStart());
				out.print(' ');
				out.print(rec.getReadLength()-(rec.getUnclippedEnd()-rec.getAlignmentEnd()));
				out.print(' ');
				out.print(rec.getReadNegativeStrandFlag()?"-":"+");
				out.print(' ');
				out.print(rec.getMappingQuality());
				out.println();
				out.println(refseq);
				out.println(readseq);
				out.println();
				
				++align_id;
				}
			out.flush();
			return 0;
			}
		catch(final Throwable err)
			{
			getLogger().error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(indexedFastaSequenceFile);
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		new Biostar170742().instanceMainWithExit(args);
		}
	}
