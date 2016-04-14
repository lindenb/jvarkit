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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.extendref;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReaderFactory.Option;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class ExtendReferenceWithReads extends AbstractExtendReferenceWithReads
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(ExtendReferenceWithReads.class);
	private List<SamReader> samReaders=new ArrayList<>();
	private enum Rescue {LEFT,CENTER,RIGHT};
	
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	

	private char consensus(Counter<Byte> count)
		{
		if(count==null || count.isEmpty()) return 'n';
		double total = (double)count.getTotal();
		if(total<  (double)this.minDepth) return 'n';
		Byte base = count.getMostFrequent();
		double nbases= (double)count.count(base);
		
		if( nbases/total >= this.callingFraction) return (char)Character.toLowerCase(base.byteValue());
		return 'n';
		}	

	/**
	 *scanRegion
	 * @param contig    Reference sequence of interest.
     * @param start     0-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
     * @param end       0-based, exclusive end of interval of interest. Zero implies end of the reference sequence.

	 */
	private Map<Integer, Counter<Byte>> scanRegion(
			final SAMSequenceRecord contig,
			final int chromStart,
			final int chromEnd,
			final Rescue rescue
			)
		{
		final Map<Integer, Counter<Byte>> pos2bases = new HashMap<>(1+chromEnd-chromStart);
		LOG.info("Scanning :"+contig.getSequenceName()+":"+chromStart+"-"+chromEnd);
		for(int side=0;side<2;++side)
			{
			if(!rescue.equals(Rescue.CENTER) && side>0)//5' or 3'
				{
				break;//already done
				}
			for(final SamReader samReader: samReaders)
				{
				final SAMSequenceDictionary dict2=samReader.getFileHeader().getSequenceDictionary();
				final SAMSequenceRecord ssr2 = dict2.getSequence(contig.getSequenceName());
				if(ssr2==null || ssr2.getSequenceLength()!=contig.getSequenceLength())
					{
					LOG.info("No contig "+contig.getSequenceName()+" with L="+contig.getSequenceLength() +" bases in "+samReader.getResourceDescription());
					continue;
					}
				int mappedPos=-1;
				
				
				switch(rescue)
					{
					case LEFT:  mappedPos = 1;break;
					case RIGHT: mappedPos = contig.getSequenceLength();break;
					case CENTER: mappedPos= (side==0? chromStart+1:chromEnd);break;
					default: throw new IllegalStateException();
					}
				final SAMRecordIterator iter =  samReader.query(
						contig.getSequenceName(),
						mappedPos,
						mappedPos,
						false
						);
				while(iter.hasNext())
					{
					final SAMRecord rec = iter.next();
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getMappingQuality()==0) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					final Cigar cigar=rec.getCigar();
					if(cigar==null) continue;
					final byte bases[]=rec.getReadBases();
					if(bases==null || bases.length==0) continue;
					int refPos1 = rec.getUnclippedStart();
					int readpos = 0;
					for(final CigarElement ce : cigar)
						{
						final CigarOperator op= ce.getOperator();
						for(int L=0;L<ce.getLength();++L)
							{
							if(op.consumesReadBases())
								{
								Counter<Byte> count= pos2bases.get(refPos1-1);
								if(count==null)
									{
									count=new Counter<>();
									pos2bases.put(refPos1-1,count);
									}
								count.incr((byte)Character.toLowerCase(bases[readpos]));
								readpos++;
								}
			
							if(op.consumesReferenceBases())
								{
								refPos1++;
								}
							}
						}
					}
				
				
				iter.close();
				}
			}
		return pos2bases;
		}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		if(super.faidx==null)
			{
			return wrapException("No REF defined option -"+OPTION_FAIDX);
			}
		final List<String> args = super.getInputFiles();
		this.samReaders.clear();
		PrintStream out=null;
		try {
			this.indexedFastaSequenceFile= new IndexedFastaSequenceFile(faidx);
			SAMSequenceDictionary dict = this.indexedFastaSequenceFile.getSequenceDictionary();
			if(dict==null)
				{
				return wrapException("Reference file is missing a sequence dictionary (use picard)");
				}
			final SamReaderFactory srf= super.createSamReaderFactory();
			srf.setOption(Option.CACHE_FILE_BASED_INDEXES, true);
			for(final String uri:IOUtils.unrollFiles(args))
				{
				LOG.info("opening BAM "+uri);
				final SamReader sr = srf.open(SamInputResource.of(uri));
					
				/* doesn't work with remote ?? */
				if(!sr.hasIndex()) 
					{
					return wrapException("file "+uri+" is not indexed");
					}
				final SAMFileHeader header= sr.getFileHeader();
				if(!header.getSortOrder().equals(SortOrder.coordinate))
					{
					return wrapException("file "+uri+" not sorted on coordinate but "+header.getSortOrder());
					}
				final SAMSequenceDictionary dict2 = header.getSequenceDictionary();
				if(dict2==null)
					{
					return wrapException("file "+uri+" needs a sequence dictionary");
					}
				
				samReaders.add(sr);
				}
			if(samReaders.isEmpty())
				{
				return wrapException("No BAM defined");
				}
			
			out = super.openFileOrStdoutAsPrintStream();
			
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(dict);
			for(final SAMSequenceRecord ssr: dict.getSequences())
				{
				final GenomicSequence genomic = new GenomicSequence(this.indexedFastaSequenceFile, ssr.getSequenceName());
				int chromStart=0;
				
				int nPrinted=0;
				out.print(">");
				out.print(ssr.getSequenceName());
				
				for(final Rescue side: Rescue.values())
					{
					switch(side) /* look 5' */
						{
						case LEFT: /* look before position 0 of chromosome */
							{
							final Map<Integer, Counter<Byte>> pos2bases = scanRegion(
									ssr,
									-1,
									-1,
									side
									);
							int newstart=0;
							for(final Integer pos:pos2bases.keySet())
								{
								if(pos>=0) continue;
								newstart=Math.min(newstart, pos);
								}
							while(newstart<0)
								{
								final Counter<Byte> count = pos2bases.get(newstart); 
								if(nPrinted%60==0) out.println();
								out.print(consensus(count));
								newstart++;
								++nPrinted;
								}
							break;
							}
						case RIGHT: /* look after position > length(chromosome) */
							{
							final Map<Integer, Counter<Byte>> pos2bases =  this.scanRegion(
									ssr,
									-1,
									-1,
									side
									);
							int newend=ssr.getSequenceLength();
							for(final Integer pos:pos2bases.keySet())
								{
								if(pos<ssr.getSequenceLength()) continue;
								newend=Math.max(newend, pos+1);
								}
							for(int i=ssr.getSequenceLength();i< newend;i++)
								{
								final Counter<Byte> count = pos2bases.get(i); 
								if(nPrinted%60==0) out.println();
								out.print(consensus(count));
								++nPrinted;
								}
							break;
							}
						case CENTER: /* 0 to chromosome-length */
							{
							chromStart=0;
							while(chromStart< genomic.length())
								{
								final char base = Character.toUpperCase(genomic.charAt(chromStart));
								if(base!='N')
									{
									progress.watch(ssr.getSequenceName(), chromStart);
									if(nPrinted%60==0) out.println();
									out.print(base);
									++chromStart;
									++nPrinted;
									continue;
									}
								
								int chromEnd=chromStart+1;
								while(chromEnd< genomic.length() && Character.toUpperCase(genomic.charAt(chromEnd))=='N')
									{
									++chromEnd;
									}
								
								if(chromEnd-chromStart< this.minLenNNNNContig)
									{
									while(chromStart<chromEnd)
										{
										progress.watch(ssr.getSequenceName(), chromStart);
										if(nPrinted%60==0) out.println();
										out.print(base);
										++chromStart;
										++nPrinted;
										}
									continue;
									}
								final Map<Integer, Counter<Byte>> pos2bases = scanRegion(
										ssr,
										chromStart,
										chromEnd,
										side
										);
								

								while(chromStart<chromEnd)
									{
									final Counter<Byte> count = pos2bases.get(chromStart); 
									if(nPrinted%60==0) out.println();
									if(count==null)
										{
										out.print('N');
										}
									else
										{
										out.print(consensus(count));
										}
									++chromStart;
									++nPrinted;
									continue;
									}
								}
							break;	
							}
						}//end switch type
					}
				out.println();
				}
			
			progress.finish();
			
			out.flush();
			out.close();
			return RETURN_OK;
			}
		catch (final Exception e) {
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			for(final SamReader r : samReaders) CloserUtil.close(r);
			samReaders.clear();
			}
		
		}
	
	public static void main(String[] args)
		{
		new ExtendReferenceWithReads().instanceMainWithExit(args);
		}

}
