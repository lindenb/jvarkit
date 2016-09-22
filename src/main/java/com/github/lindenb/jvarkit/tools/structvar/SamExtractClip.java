/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar;

import java.util.Collection;
import java.util.List;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class SamExtractClip extends AbstractSamExtractClip
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SamExtractClip.class);

	@Override
	public Collection<Throwable> call() throws Exception {
		SamReader r=null;
		BasicFastqWriter out=null;
		final List<String> args = super.getInputFiles();
		try
				{
				if(getOutputFile()!=null)
					{
					LOG.info("writing to "+getOutputFile());
					out=new BasicFastqWriter(getOutputFile());
					}
				else
					{
					LOG.info("writing to stdout");
					out=new BasicFastqWriter(stdout());
					}
				if(args.isEmpty())
					{
					LOG.info("Reading from stdin");
					r= createSamReaderFactory().open(SamInputResource.of(stdin()));
					run(r,out);
					r.close();
					}
				else 
					{
					for(final String filename:args)
						{
						LOG.info("Reading from "+filename);
						r= createSamReaderFactory().open(SamInputResource.of(filename));
						run(r,out);
						r.close();r=null;
						}
					}
				out.flush();
				return RETURN_OK;
				}
			catch(final Exception err)
				{
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(r);
				CloserUtil.close(out);
				}
		}
		
		private void run(final SamReader r,final FastqWriter out)
			{
			int startend[]=new int[2];
			final SAMFileHeader header=r.getFileHeader();
			//w=swf.make(header, System.out);
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header);
			SAMRecordIterator it= r.iterator();
			while(it.hasNext())
				{
				final SAMRecord rec=progress.watch(it.next());
				
				if(rec.getReadUnmappedFlag()) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				
				final Cigar cigar=rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				
				String suffix="";
				if(rec.getReadPairedFlag())
					{
					suffix=(rec.getFirstOfPairFlag()?"/1":"/2");
					}
				
			
				
				
				startend[0]=0;
				startend[1]=rec.getReadLength();
				boolean found=false;
				for(int side=0;side<2;++side)
					{
					final CigarElement ce=cigar.getCigarElement(side==0?0:cigar.numCigarElements()-1);
					if(!ce.getOperator().equals(CigarOperator.S)) continue;
					if(ce.getLength() < min_clip_length) continue;
					
					found=true;
					final String clippedSeq;
					final String clippedQual;
					if(side==0)
						{
						startend[0]=ce.getLength();
						clippedSeq= rec.getReadString().substring(0, startend[0]);
						clippedQual=rec.getBaseQualityString().substring(0, startend[0]);
						}
					else
						{
						startend[1]=rec.getReadLength()-ce.getLength();
						clippedSeq= rec.getReadString().substring(startend[1]);
						clippedQual=rec.getBaseQualityString().substring(startend[1]);
						}
					
					out.write(new FastqRecord(
							rec.getReadName()+suffix+";"+side+";"+rec.getReferenceName()+";"+rec.getAlignmentStart()+";"+rec.getFlags()+";"+rec.getCigarString()+";"+(side==0?"5'":"3'"),
							clippedSeq,
							"",
							clippedQual
							));
					}
				if(!found) continue;
				
				String bases=rec.getReadString();
				String qual=rec.getBaseQualityString();
				if( rec.getReadNegativeStrandFlag())
					{
					bases=AcidNucleics.reverseComplement(bases);
					qual=new StringBuilder(qual).reverse().toString();
					}
				
				if(print_original_read)
					{
					out.write(new FastqRecord(
							rec.getReadName()+suffix,
							bases,
							"",
							qual
							));
					}
				
				if(print_clipped_read)
					{
					out.write(new FastqRecord(
							rec.getReadName()+suffix+":clipped",
							bases.substring(startend[0], startend[1]),
							"",
							qual.substring(startend[0], startend[1])
							));
					}
				}
			
			it.close();
			progress.finish();
			}

		
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamExtractClip().instanceMainWithExit(args);
		}
	}
