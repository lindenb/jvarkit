
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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.bamstats01;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.SequenceUtil;

/**

BEGIN_DOC


### History

* Dec 2013 Added PROPER_PAIR_HMQ for @SolenaLS
* Dec 2013 Added X and Y for @SolenaLS


### Output

See also: http://picard.sourceforge.net/explain-flags.html


#### Counts


* TOTAL : total number of reads (not PAIRS of reads)
* PAIRED: total number of reads paired (should be equals to ALL for a paired-end assay)
* UNMAPPED : count unmapped reads 
* MAPPED  : count mapped reads
* PROPER_PAIR  : count reads in proper pair (forward+reverse, same chromosome, distance is OK)
* PROPER_PAIR_HMQ  : proper pairs with mapping quality >= user qual
* PLUS : reads on plus strand
* MINUS : reads on minus strand
* PRIMARY_ALIGNMENT : alignment flagged as primary alignment (not alternative position)
* FAIL_MAPPING_QUALITY : MAQ < user qual
* DUPLICATE : the flag 'duplicate' was set
* FAIL_VENDOR_QUALITY : the flag "read fails platform/vendor quality checks" was set
* OK_FOR_PE_CALLING : reads ok for Paired-end mapping ( properly paired, not dup, not fails_vendor_qual,  not fails_mapping_qual, primary align )
* X and Y : number of reads mapping the chromosomes X/chrX and Y/chrY


#### Categories


* ALL: all reads
* IN_TARGET: reads overlapping user's BED (if provided)
* OFF_TARGET: reads with no overlap with user's BED (if provided)



### Example


```
$  java -jar dist/bamstats01.jar \
		-B capture.bed my.bam \
		

(...)
#Filename	Sample	ALL_TOTAL	ALL_PAIRED	ALL_UNMAPPED	ALL_MAPPED	ALL_PROPER_PAIR	ALL_PLUS_STRAND	ALL_MINUS_STRAND	ALL_PRIMARY_ALIGNMENT	ALL_FAIL_MAPPING_QUALITY	ALL_DUPLICATE	ALL_FAIL_VENDOR_QUALITY	IN_TARGET_TOTAL	IN_TARGET_PAIRED	IN_TARGET_UNMAPPED	IN_TARGET_MAPPED	IN_TARGET_PROPER_PAIR	IN_TARGET_PLUS_STRAND	IN_TARGET_MINUS_STRAND	IN_TARGET_PRIMARY_ALIGNMENT	IN_TARGET_FAIL_MAPPING_QUALITY	IN_TARGET_DUPLICATE	IN_TARGET_FAIL_VENDOR_QUALITY	OFF_TARGET_TOTAL	OFF_TARGET_PAIRED	OFF_TARGET_UNMAPPED	OFF_TARGET_MAPPED	OFF_TARGET_PROPER_PAIR	OFF_TARGET_PLUS_STRAND	OFF_TARGET_MINUS_STRAND	OFF_TARGET_PRIMARY_ALIGNMENT	OFF_TARGET_FAIL_MAPPING_QUALITY	OFF_TARGET_DUPLICATE	OFF_TARGET_FAIL_VENDOR_QUALITY
my.bam	Sample	1617984	1617984	3966	1614018	1407862	806964	807054	1614018	56980	0	0	1293922	1293922	0	1293922	1133808	644741	649181	1293922	14087	0	0	320096	320096	0	320096	274054	162223	157873	320096	42893	0	0
(...)

```


END_DOC
*/


@Program(name="samstats01",
	description="Statistics about the reads in a BAM.",
	keywords= {"sam","bam"},
	modificationDate="20190506"
	)
public class BamStats01
	extends Launcher
	{

	private static final Logger LOG = Logger.build(BamStats01.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	@Parameter(names={"-B","--bed"},description="capture bed file. Optional")
	private Path bedFile = null;

	@Parameter(names={"-q","--qual"},description="min mapping quality")
	private double minMappingQuality = 30.0 ;
	
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition groupBy=SAMRecordPartition.sample;

	

	private PrintStream out=System.out;
	//private File bedFile=null;
	//private double minMappingQuality=30.0;
	private int chrX_index=-1;
	private int chrY_index=-1;
	private SAMSequenceDictionary samSequenceDictionary=null;
	private IntervalTreeMap<Boolean> intervalTreeMap=null;
    
    private  class Histogram2
    	{
    	long counts[];
    	
    	Histogram2()
    		{
    		counts=new long[Category.values().length];
    		Arrays.fill(counts, 0L);
    		}
    	void increment(Category cat)
    		{
    		counts[cat.ordinal()]++;
    		}
    	void watch(SAMRecord rec)
    		{
    		boolean ok_pe_alignment=true;
			this.increment(Category.TOTAL);
			
			if(rec.getReadPairedFlag())
				{
				this.increment(Category.PAIRED);
				}
			else
				{
				ok_pe_alignment=false;
				}
			
			if(rec.getReadFailsVendorQualityCheckFlag())
				{
				ok_pe_alignment=false;
				this.increment(Category.FAIL_VENDOR_QUALITY);
				}
			
			/* count reads = unmapped+ primary align */
			if(rec.getReadUnmappedFlag() ||
				!rec.isSecondaryOrSupplementary())
				{
				this.increment(Category.UNMAPPED_PLUS_PRIMARY);
				}
			
			
			if(rec.getReadUnmappedFlag())
				{
				ok_pe_alignment=false;
				this.increment(Category.UNMAPPED);
				return;
				}
			/* mapped below ************************************/
			
			this.increment(Category.MAPPED);
			
			if(rec.getReadPairedFlag() && rec.getProperPairFlag())
				{
				this.increment(Category.PROPER_PAIR);
				if(rec.getMappingQuality()>=minMappingQuality &&
						rec.getMappingQuality()!=0 &&
						rec.getMappingQuality()!=255)
					{
					this.increment(Category.PROPER_PAIR_HMQ);
					}
				}
			else
				{
				ok_pe_alignment=false;	
				}		
			
			if(rec.getMappingQuality()==0 ||
				rec.getMappingQuality()==255)
				{
				this.increment(Category.BAD_MAPQ);
				ok_pe_alignment=false;	
				}
			
			if(rec.getReadNegativeStrandFlag())
				{
				this.increment(Category.MINUS_STRAND);
				}
			else
				{
				this.increment(Category.PLUS_STRAND);
				}
			
			
			if(!rec.isSecondaryAlignment())
				{
				this.increment(Category.PRIMARY_ALIGNMENT);
				}
			else
				{
				ok_pe_alignment=false;
				}
			
			if( !rec.isSecondaryOrSupplementary() && // bwa mem -P
				rec.getDuplicateReadFlag())
				{
				this.increment(Category.DUPLICATE);
				ok_pe_alignment=false;
				}
			
			if(rec.getMappingQuality()<minMappingQuality)
				{
				this.increment(Category.FAIL_MAPPING_QUALITY);
				ok_pe_alignment=false;
				}
			
			if(rec.getSupplementaryAlignmentFlag())
				{
				this.increment(Category.SUPPLEMENTARY_ALIGNMENT);
				ok_pe_alignment=false;
				}
			
			if(ok_pe_alignment)
				{
				this.increment(Category.OK_FOR_PE_CALLING);
				if(chrX_index==rec.getReferenceIndex())
					{
					this.increment(Category.X);
					}
				else if(chrY_index==rec.getReferenceIndex())
					{
					this.increment(Category.Y);
					}

				}
			
			
    		}
    	
    	}
    
    private  class Histogram
    	{
    	final Histogram2 histograms[]=new Histogram2[]{
    			new Histogram2(),//ALL
    			new Histogram2(),//in capture
    			new Histogram2()//off capture
    			};
    	}

    
    private enum Category2
    	{
    	ALL,IN_TARGET,OFF_TARGET
    	}
	
	
	private enum Category
		{
		TOTAL,
		UNMAPPED_PLUS_PRIMARY,
		PAIRED,
		UNMAPPED,
		MAPPED,
		PROPER_PAIR,
		PROPER_PAIR_HMQ,
		PLUS_STRAND,
		MINUS_STRAND,
		PRIMARY_ALIGNMENT,
		FAIL_MAPPING_QUALITY,
		DUPLICATE,
		FAIL_VENDOR_QUALITY,
		OK_FOR_PE_CALLING,
		SUPPLEMENTARY_ALIGNMENT,//new in samSpec 2013-12
		X,Y,
		BAD_MAPQ
		;
		};
	
	private BamStats01()
		{
		
		}
		
	private void run(final String filename,SamReader samFileReader) throws IOException
		{	
		Map<String,Histogram> sample2hist=new HashMap<String, BamStats01.Histogram>();
			
			SAMSequenceDictionary currDict=samFileReader.getFileHeader().getSequenceDictionary();

			if(samSequenceDictionary==null)
				{
				samSequenceDictionary=currDict;
				}
			
			if(this.bedFile!=null )
				{
				if(!SequenceUtil.areSequenceDictionariesEqual(currDict, samSequenceDictionary))
					{
					samFileReader.close();
					throw new IOException("incompatible sequence dictionaries."+filename);
					}
					
				
				if(intervalTreeMap==null)
					{
					final BedLineCodec bedCodec=new BedLineCodec();
					intervalTreeMap=new IntervalTreeMap<>();
					LOG.info("opening "+this.bedFile);
					String line;
					final BufferedReader bedIn=IOUtils.openPathForBufferedReading(bedFile);
					while((line=bedIn.readLine())!=null)
						{
						final BedLine bedLine = bedCodec.decode(line);
						if(bedLine==null) continue;
						int seqIndex=currDict.getSequenceIndex(bedLine.getContig());
						if(seqIndex==-1)
							{
							throw new JvarkitException.ContigNotFoundInDictionary(bedLine.getContig(),currDict);
							}
						intervalTreeMap.put(bedLine.toInterval(),Boolean.TRUE);
						}
					bedIn.close();
					LOG.info("done reading "+this.bedFile);
					}
				}
			this.chrX_index=-1;
			this.chrY_index=-1;
			
			
			for(final SAMSequenceRecord rec:currDict.getSequences())
				{
				final String chromName=rec.getSequenceName().toLowerCase();
				if(chromName.equals("x") || chromName.equals("chrx"))
					{
					this.chrX_index=rec.getSequenceIndex();
					}
				else if(chromName.equals("y") || chromName.equals("chry"))
					{
					this.chrY_index=rec.getSequenceIndex();
					}
				}
			
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(currDict);
			final SAMRecordIterator iter=samFileReader.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec=progess.watch(iter.next());
				
				String sampleName = groupBy.getPartion(rec);
				if(sampleName==null || sampleName.isEmpty()) sampleName="undefined";
				
				Histogram hist=sample2hist.get(sampleName);
				if(hist==null)
					{
					hist=new Histogram();
					sample2hist.put(sampleName, hist);
					}
				
				hist.histograms[Category2.ALL.ordinal()].watch(rec);
				
				
				
				if(intervalTreeMap==null) continue;
				if(rec.getReadUnmappedFlag())
					{
					continue;
					}
				
			
				
				if(!intervalTreeMap.containsOverlapping(new Interval(
							rec.getReferenceName(),
							rec.getAlignmentStart(),
							rec.getAlignmentEnd()
							)))
					{
					hist.histograms[Category2.OFF_TARGET.ordinal()].watch(rec);
					}		
				else
					{
					hist.histograms[Category2.IN_TARGET.ordinal()].watch(rec);
					}
				}
			progess.finish();
			samFileReader.close();
			samFileReader=null;
		
			for(final String sampleName: sample2hist.keySet())
				{
				final Histogram hist=sample2hist.get(sampleName);
				out.print(filename+"\t"+sampleName);
				
				for(final Category2 cat2: Category2.values())
					{
					for(final Category cat1: Category.values())//je je suis libertineuuh, je suis une cat1
						{
						out.print("\t");
						out.print(hist.histograms[cat2.ordinal()].counts[cat1.ordinal()]);
						}
					if(intervalTreeMap==null) break;
					}
				out.println();
				}
		}
	@Override
	public int doWork(final List<String> inputs) {
		final List<String> args= new ArrayList<>(IOUtils.unrollFiles(inputs)); 
		try {
			this.out = super.openPathOrStdoutAsPrintStream(this.outputFile);
			
			out.print("#Filename\tSample");
			for(Category2 cat2: Category2.values())
				{
				for(Category cat1: Category.values())
					{
					out.print("\t"+cat2+"_"+cat1);//je je suis libertineuuh, je suis une cat1
					}
				if(bedFile==null) break;
				}
			out.println();
			
			final SamReaderFactory srf= super.createSamReaderFactory();
			
			if(args.isEmpty())
				{
				final SamReader r= srf.open(SamInputResource.of(stdin()));
				run("stdin",r);
				r.close();
				}
			else
				{
				for(final String filename:args)
					{
					try(final SamReader sfr=srf.open(Paths.get(filename)))
						{
						run(filename,sfr);
						}
					}
				}
			out.flush();			
			this.out.flush();
			this.out.close();
			return RETURN_OK;
		} catch (final Exception e) {
			LOG.error(e);
			return -1;
		} finally
			{
			CloserUtil.close(out);
			}
		}
	
	
	public static void main(final String[] args)
		{
		new BamStats01().instanceMainWithExit(args);
		}
	}
