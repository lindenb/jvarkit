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
package com.github.lindenb.jvarkit.tools.vcfannobam;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SequenceUtil;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## Example

```bash
$  java -jar dist/vcfannobam.jar \
		-BAM input.bam\
		-BED capture.bed \
		input.vcf.gz

(...)
##INFO=<ID=CAPTURE,Number=1,Type=String,Description="Capture stats: Format is (start|end|mean|min|max|length|not_covered|percent_covered) ">
(...)
2	16100665	.	A	T	13261.77	.	CAPTURE=16100619|16100715|1331.96|1026.0|1773.0|97|0|100
2	178395141	.	T	A	1940.77	.	CAPTURE=178394991|178395199|193.11|100.0|276.0|209|0|100
(...)
```

END_DOC

 */
@Program(name="vcfannobam",
	deprecatedMsg="useless: use DP/DP4 in the Genotypes, or use GATK variant annotator",
	description="Annotate a VCF with the Coverage statistics of a BAM file+  BED file of capture. ",
	keywords= {"bam","sam","depth","vcf"})
public class VCFAnnoBam extends Launcher {
	private static final Logger LOG = Logger.build(VCFAnnoBam.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
    @Parameter(names= {"-BED","--bed"}, description="BED File Exome capture.",required=true)
	private File BEDILE=null;
	@Parameter(names= {"-BAM","--bam"}, description="Indexed BAM File. One file ending with the '.list' suffix will be interpreted as a text file containing the path to the bams.")		
	private List<File> BAMFILE=null;
	@Parameter(names={"-filter","--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names= {"-MIN_COV","--coverage"}, description="min coverage to say the position is not covered")		
	private int MIN_COVERAGE=0;
	@Parameter(names= {"-tag","--tag"}, description="VCF info tag")		
	private String capture_tag="CAPTURE";
    
    private class Rgn
    	{
    	Interval interval;
 			
    		int count_no_coverage=0;
    		double mean=0;
        	double min=0.0;
        	double max=-1;
        	int percent_covered=0;
        	
       
    	
    	boolean processed=false;
    	
    	@Override
    	public String toString()
    		{
    		final StringBuilder b=new StringBuilder();
    		b.append(interval.getStart());
    		b.append('|');
    		b.append(interval.getEnd());
    		b.append('|');
    		b.append(String.format("%.2f",mean));
    		b.append('|');
    		b.append(min);
    		b.append('|');
    		b.append(max);
    		b.append('|');
    		b.append((interval.getEnd()-interval.getStart()+1));
    		b.append('|');
    		b.append(count_no_coverage);
    		b.append('|');
    		b.append(percent_covered);
    		return b.toString();
    		}

    	}

    private void process(Rgn rgn,List<SamReader> samReaders)
		{
    	rgn.processed=true;
		int chromStart1= rgn.interval.getStart();
		int chromEnd1=  rgn.interval.getEnd();
		
		
		int counts[]=new int[chromEnd1-chromStart1+1];
		if(counts.length==0) return;
		Arrays.fill(counts, 0);
		
		for(SamReader samReader:samReaders)
			{
			/**
			 *     start - 1-based, inclusive start of interval of interest. Zero implies start of the reference sequence.
			*	   end - 1-based, inclusive end of interval of interest. Zero implies end of the reference sequence. 
			 */
			SAMRecordIterator r=samReader.queryOverlapping(rgn.interval.getContig(), chromStart1, chromEnd1);
			while(r.hasNext())
				{
				SAMRecord rec=r.next();
				if(rec.getReadUnmappedFlag()) continue;
				if(this.filter.filterOut(rec)) continue;
				if(!rec.getReferenceName().equals(rgn.interval.getContig())) continue;
			
				Cigar cigar=rec.getCigar();
				if(cigar==null) continue;
	    		int refpos1=rec.getAlignmentStart();
	    		for(CigarElement ce:cigar.getCigarElements())
	    			{
					switch(ce.getOperator())
						{
						case H:break;
						case S:break;
						case I:break;
						case P:break;
						case N:// reference skip
						case D://deletion in reference
							{
	    					refpos1+=ce.getLength();
							break;
							}
						case M:
						case EQ:
						case X:
							{
							for(int i=0;i< ce.getLength() && refpos1<= chromEnd1;++i)
	    		    			{
								if(refpos1>= chromStart1 && refpos1<=chromEnd1)
									{
									counts[refpos1-chromStart1]++;
									}
	    						refpos1++;
			    				}
							break;
							}
						default: throw new IllegalStateException(
								"Doesn't know how to handle cigar operator:"+ce.getOperator()+
								" cigar:"+cigar
								);
	
						}
	    				
	    			}
				}
			
			r.close();
			}
		
			Arrays.sort(counts);
			
			for(int cov:counts)
				{
				if(cov<=MIN_COVERAGE) rgn.count_no_coverage++;
				rgn.mean+=cov;
				}
			rgn.mean/=counts.length;
			rgn.min=counts[0];
			rgn.max=counts[counts.length-1];
			rgn.percent_covered=(int)(((counts.length-rgn.count_no_coverage)/(double)counts.length)*100.0);
			rgn.processed=true;
			
			}
	
    @Override
    protected int doVcfToVcf(final String inputName, final VCFIterator r,final  VariantContextWriter w) {
		BufferedReader bedIn=null;
		final List<SamReader> samReaders=new ArrayList<SamReader>();
		final IntervalTreeMap<Rgn> capture=new IntervalTreeMap<Rgn>();
		try
			{
			SAMFileHeader firstHeader=null;
			
			if(this.BAMFILE.size()==1 && this.BAMFILE.get(0).getName().endsWith(".list"))
				{
				final File first = this.BAMFILE.get(0);
				this.BAMFILE.clear();
				this.BAMFILE.addAll(IOUtils.unrollFile(first));
				}
			
			for(final  File samFile:new HashSet<File>(BAMFILE))
				{
				LOG.info("open bam "+samFile);
				final SamReader samReader = super.openSamReader(samFile.getPath());
				final SAMFileHeader samHeader=samReader.getFileHeader();
				samReaders.add(samReader);
				if(firstHeader==null)
					{
					firstHeader=samHeader;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(
						firstHeader.getSequenceDictionary(),
						samHeader.getSequenceDictionary())
						)
					{
					throw new JvarkitException.DictionariesAreNotTheSame(
							firstHeader.getSequenceDictionary(),
							samHeader.getSequenceDictionary());
					}
				}		
			IntervalList intervalList=new IntervalList(firstHeader);
			LOG.info("read bed "+BEDILE);
			bedIn=IOUtils.openFileForBufferedReading(BEDILE);
			String line;
			final BedLineCodec bedCodec=new BedLineCodec();
			while((line=bedIn.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				final BedLine bed=bedCodec.decode(line);
				if(bed==null) continue;
				
				if(firstHeader.getSequenceDictionary().getSequence(bed.getContig())==null)
					{
					LOG.error("error in BED +"+BEDILE+" : "+line+" chromosome is not in sequence dict of "+BAMFILE);
					continue;
					}
				intervalList.add(bed.toInterval());
				}
			bedIn.close();
			bedIn=null;
			intervalList=intervalList.sorted();
			for(final Interval interval:intervalList.uniqued())
				{				
				final Rgn rgn=new Rgn();
				rgn.interval=interval;
				capture.put(rgn.interval, rgn);
				}
			intervalList=null;
			
			VCFHeader header=r.getHeader();
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			h2.addMetaDataLine(new VCFInfoHeaderLine(
					this.capture_tag, 1,
					VCFHeaderLineType.String,
					"Capture stats: Format is (start|end|mean|min|max|length|not_covered|percent_covered) BAM files: "+BAMFILE+" CAPTURE:"+BEDILE));
			w.writeHeader(h2);
			
			
			while(r.hasNext())
				{
				final VariantContext ctx=r.next();
				Interval interval=new Interval(ctx.getContig(), ctx.getStart(), ctx.getEnd());
				Collection<Rgn> rgns=capture.getOverlapping(interval);
				Iterator<Rgn> it=rgns.iterator();
				if(!it.hasNext())
					{
					w.add(ctx);
					continue;
					}
				final Rgn rgn=it.next();
				if(!rgn.processed)
					{
					//LOG.info("processing "+rgn.interval);
					process(rgn,samReaders);
					}
				final VariantContextBuilder b=new VariantContextBuilder(ctx);
				b.attribute(this.capture_tag, rgn.toString());
				w.add(b.make());
				}
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			for(final SamReader samReader:samReaders) CloserUtil.close(samReader);
			}
		}

    @Override
    public int doWork(final List<String> args) {
    	return doVcfToVcf(args, outputFile);
    	}
    
	/**
	 * @param args
	 */
	public static void main(final String[] args) {
		new VCFAnnoBam().instanceMainWithExit(args);
	}

}
