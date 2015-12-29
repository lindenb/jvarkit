/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
import java.util.Collection;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;


public class Biostar170742 extends AbstractBiostar170742
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(Biostar170742.class);

	
	
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(super.faidx==null)
			{
			return wrapException("Reference sequence was not defined");
			}
		PrintStream out=null;
		SamReader sfr=null;
		SAMRecordIterator iter=null;
		GenomicSequence genomicSequence = null;
		IndexedFastaSequenceFile indexedFastaSequenceFile= null;
		try
			{
			indexedFastaSequenceFile = new IndexedFastaSequenceFile(super.faidx);
			long align_id=0;
			sfr = openSamReader(inputName);
			out =  super.openFileOrStdoutAsPrintStream();
			final StringBuilder refseq = new StringBuilder();
			final StringBuilder readseq = new StringBuilder();
			
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(sfr.getFileHeader());
			iter=sfr.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec= progress.watch(iter.next());
				if(rec.getReadUnmappedFlag()) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar==null) continue;
				final byte readbases[] = rec.getReadBases();
				if(readbases==null) continue;

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
			progress.finish();
			iter.close();
			out.flush();
			LOG.info("done");
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(iter);
			CloserUtil.close(sfr);
			CloserUtil.close(indexedFastaSequenceFile);
			}
		}
	
	public static void main(String[] args)throws Exception
		{
		new Biostar170742().instanceMain(args);
		}
	}
