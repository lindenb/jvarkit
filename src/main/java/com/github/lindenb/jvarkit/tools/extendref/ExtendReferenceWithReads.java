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
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class ExtendReferenceWithReads extends AbstractCommandLineProgram
	{
	private int minLenNNNNContig=100;
	private double callingFraction=0.8;
	private int minDepth=1;
	private List<SamReader> samReaders=new ArrayList<>();
	private enum Rescue {LEFT,CENTER,RIGHT};
	
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	@Override
	protected String getOnlineDocUrl()
		{
		return DEFAULT_WIKI_PREFIX+"ExtendReferenceWithReads";
		}
	@Override
	public String getProgramDescription() {
		return "Extending ends of sequences with the help of reads https://www.biostars.org/p/148089/";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -R (ref) "+getMessageBundle("reference.faidx")+" Required");
		out.println(" -N (int) consider only gaps with size>=N default:"+this.minLenNNNNContig);
		out.println(" -f (0.0<float<=1.0) new base must have fraction greater than this number :"+this.callingFraction);
		out.println(" -d (int) min depth. default:"+this.minDepth);
		super.printOptions(out);
		}

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
			SAMSequenceRecord contig,
			int chromStart,
			int chromEnd,
			Rescue rescue
			)
		{
		Map<Integer, Counter<Byte>> pos2bases = new HashMap<>(1+chromEnd-chromStart);
		info("Scanning :"+contig.getSequenceName()+":"+chromStart+"-"+chromEnd);
		for(int side=0;side<2;++side)
			{
			if(!rescue.equals(Rescue.CENTER) && side>0)//5' or 3'
				{
				break;//already done
				}
			for(SamReader samReader: samReaders)
				{
				SAMSequenceDictionary dict2=samReader.getFileHeader().getSequenceDictionary();
				SAMSequenceRecord ssr2 = dict2.getSequence(contig.getSequenceName());
				if(ssr2==null || ssr2.getSequenceLength()!=contig.getSequenceLength())
					{
					warning("No contig "+contig.getSequenceName()+" with L="+contig.getSequenceLength() +" bases in "+samReader.getResourceDescription());
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
				SAMRecordIterator iter =  samReader.query(
						contig.getSequenceName(),
						mappedPos,
						mappedPos,
						false
						);
				while(iter.hasNext())
					{
					SAMRecord rec = iter.next();
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.getMappingQuality()==0) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					Cigar cigar=rec.getCigar();
					if(cigar==null) continue;
					byte bases[]=rec.getReadBases();
					if(bases==null || bases.length==0) continue;
					int refPos1 = rec.getUnclippedStart();
					int readpos = 0;
					for(CigarElement ce: cigar.getCigarElements())
						{
						CigarOperator op= ce.getOperator();
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
	public int doWork(String[] args)
		{
		File faidx=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"R:N:f:d:"))!=-1)
			{
			switch(c)
				{
				case 'R': faidx = new File(opt.getOptArg());break;
				case 'N': this.minLenNNNNContig = Math.max(1, Integer.parseInt(opt.getOptArg()));break;
				case 'f': this.callingFraction = Math.min(1.0,Math.max(0.0, Double.parseDouble(opt.getOptArg())));break;
				case 'd': this.minDepth = Math.max(1, Integer.parseInt(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(faidx==null)
			{
			error("No REF defined");
			return -1;
			}
		this.samReaders.clear();
		PrintStream out=System.out;
		try {
			this.indexedFastaSequenceFile= new IndexedFastaSequenceFile(faidx);
			SAMSequenceDictionary dict = this.indexedFastaSequenceFile.getSequenceDictionary();
			if(dict==null)
				{
				error("Reference file is missing a sequence dictionary (use picard)");
				return -1;
				}
			SamReaderFactory srf=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			srf.setOption(Option.CACHE_FILE_BASED_INDEXES, true);
			for(String uri:IOUtils.unrollFiles(Arrays.asList(Arrays.copyOfRange(args, opt.getOptInd(), args.length))))
				{
				info("opening BAM "+uri);
				SamReader sr = null;
				
				sr = srf.open(SamInputResource.of(uri));
					
				
				
				/* doesn't work with remote ?? */
				if(!sr.hasIndex()) 
					{
					error("file "+uri+" is not indexed");
					return -1;
					}
				SAMFileHeader header= sr.getFileHeader();
				if(!header.getSortOrder().equals(SortOrder.coordinate))
					{
					error("file "+uri+" not sorted on coordinate but "+header.getSortOrder());
					return -1;
					}
				SAMSequenceDictionary dict2 = header.getSequenceDictionary();
				if(dict2==null)
					{
					error("file "+uri+" needs a sequence dictionary");
					return -1;
					}
				
				samReaders.add(sr);
				}
			if(samReaders.isEmpty())
				{
				error("No BAM defined");
				return -1;
				}
			SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(dict);
			for(SAMSequenceRecord ssr: dict.getSequences())
				{
				GenomicSequence genomic=new GenomicSequence(this.indexedFastaSequenceFile, ssr.getSequenceName());
				int chromStart=0;
				
				int nPrinted=0;
				out.print(">");
				out.print(ssr.getSequenceName());
				
				for(Rescue side: Rescue.values())
					{
					switch(side) /* look 5' */
						{
						case LEFT: /* look before position 0 of chromosome */
							{
							Map<Integer, Counter<Byte>> pos2bases = scanRegion(
									ssr,
									-1,
									-1,
									side
									);
							int newstart=0;
							for(Integer pos:pos2bases.keySet())
								{
								if(pos>=0) continue;
								newstart=Math.min(newstart, pos);
								}
							while(newstart<0)
								{
								Counter<Byte> count = pos2bases.get(newstart); 
								if(nPrinted%60==0) out.println();
								out.print(consensus(count));
								newstart++;
								++nPrinted;
								}
							break;
							}
						case RIGHT: /* look before position > length(chromosome) */
							{
							Map<Integer, Counter<Byte>> pos2bases =  scanRegion(
									ssr,
									-1,
									-1,
									side
									);
							int newend=ssr.getSequenceLength();
							for(Integer pos:pos2bases.keySet())
								{
								if(pos<ssr.getSequenceLength()) continue;
								newend=Math.max(newend, pos+1);
								}
							for(int i=ssr.getSequenceLength();i< newend;i++)
								{
								Counter<Byte> count = pos2bases.get(i); 
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
								char base = Character.toUpperCase(genomic.charAt(chromStart));
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
								Map<Integer, Counter<Byte>> pos2bases = scanRegion(
										ssr,
										chromStart,
										chromEnd,
										side
										);
								

								while(chromStart<chromEnd)
									{
									Counter<Byte> count = pos2bases.get(chromStart); 
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
			return 0;
			}
		catch (Exception e) {
			error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			for(SamReader r : samReaders) CloserUtil.close(r);
			samReaders.clear();
			}
		
		}
	
	public static void main(String[] args)
		{
		new ExtendReferenceWithReads().instanceMainWithExit(args);
		}

}
