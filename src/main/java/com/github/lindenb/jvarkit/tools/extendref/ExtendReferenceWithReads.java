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
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamRecordFilterFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

/**
 BEGIN_DOC
 
 ## Example

```
$  java   -jar dist/extendrefwithreads.jar \
     -R human_g1k_v37.fasta -f 0.3 \
     f1.bam f2.bam f3.bam 2> /dev/null |\
  cat -n | grep -E '(>|[atgc])' 

     1	>1
   167	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNncgattaccctaacgctcac
   168	cctaaccctcnccctntnccnncnncccnncttcttccgaTAACCCTAACCCTAACCCTA
  3791	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNatt
  3792	tatgcNctttntgctgtGATTCATGGCTGAAATCGTGTTTGACCAGCTATGTGTGTCTCT
  8691	NNNNNNNNNNNNNNNNNNNNNNNNctagGATCCTTGAAGCGCCCCCAAGGGCATCTTCTC
 64089	TGGTGAGGGAAATTAGAACCACGACAATTTGGGAACTTAGCTTCTGCCctgctccNNNNN
 66589	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNgagtAGCTGAGACTAC
 
 ```

 END_DOC
 */
@Program(name="extendrefwithreads",
	description="Extending ends of sequences with the help of reads",
	biostars=148089,
	keywords={"read","fastq","reference","sam","bam"},
	modificationDate="20190417"
	)
public class ExtendReferenceWithReads extends Launcher
	{
	private static final Logger LOG = Logger.build(ExtendReferenceWithReads.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx=null;
	@Parameter(names={"-f","--callingfraction"},description="(0.0<float<=1.0) new base must have fraction greater than this number")
	private double callingFraction=0.8;
	@Parameter(names={"-d","--mindepth"},description="min depth")
	private int minDepth=1;
	@Parameter(names={"-N","--mincontig"},description="onsider only gaps in reference with size&gt;=N")
	private int minLenNNNNContig=100;
	@Parameter(names={"-filter","--filter"},description=SamRecordFilterFactory.FILTER_DESCRIPTION,converter=SamRecordFilterFactory.class,splitter=NoSplitter.class)
	private SamRecordFilter filter  = SamRecordFilterFactory.getDefault();

	
	
	private List<SamReader> samReaders=new ArrayList<>();
	private enum Rescue {LEFT,CENTER,RIGHT};
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	

	private char consensus(final Counter<Byte> count)
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
					if(this.filter.filterOut(rec)) continue;
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
	public int doWork(List<String> args) {
		if(this.faidx==null)
			{
			LOG.error("No REF defined");
			return -1;
			}
		this.samReaders.clear();
		PrintStream out=null;
		try {
			this.indexedFastaSequenceFile= new IndexedFastaSequenceFile(faidx);
			SAMSequenceDictionary dict = this.indexedFastaSequenceFile.getSequenceDictionary();
			if(dict==null)
				{
				LOG.error("Reference file is missing a sequence dictionary (use picard)");
				return -1;
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
					LOG.error("file "+uri+" is not indexed");
					return -1;
					}
				final SAMFileHeader header= sr.getFileHeader();
				if(!header.getSortOrder().equals(SortOrder.coordinate))
					{
					LOG.error("file "+uri+" not sorted on coordinate but "+header.getSortOrder());
					return -1;
					}
				final SAMSequenceDictionary dict2 = header.getSequenceDictionary();
				if(dict2==null)
					{
					LOG.error("file "+uri+" needs a sequence dictionary");
					return -1;
					}
				
				samReaders.add(sr);
				}
			if(samReaders.isEmpty())
				{
				LOG.error("No BAM defined");
				return -1;
				}
			
			out = super.openFileOrStdoutAsPrintStream(this.outputFile);
			
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
			LOG.error(e);
			return -1;
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
