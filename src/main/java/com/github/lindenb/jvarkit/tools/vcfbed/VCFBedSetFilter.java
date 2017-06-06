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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.File;
import java.util.List;

import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**
 * 
 * VCFBedSetFilter
 *
 */
/**

BEGIN_DOC




### Examples

```
$java -jar dist/vcfbedsetfilter.jar -f MYFILTER - -B in.bed in.vcf 

```

END_DOC
*/



@Program(name="vcfbedsetfilter",
	description="Set FILTER for VCF if it doesn't intersects with BED.",
	deprecatedMsg="use GATK FilterVariants",
	keywords={"vcf","bed","filter"}
		)
public class VCFBedSetFilter extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFBedSetFilter.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-f","--filter"},description="FILTER name")
	private String filterName = "VCFBED";

	@Parameter(names={"-i","--inverse"},description="inverse selection")
	private boolean inverse = false;

	@Parameter(names={"-B","--bed"},description="Tribble or Tabix bed file")
	private File tabixFile = null;

	@Parameter(names={"-m","--map"},description="unindexed bed file, will be loaded in memory (faster than tribble/tabix but memory consumming)")
	private File treeMapFile = null;

	@Parameter(names={"-d","--discard"},description="Discard filtered variants")
	private boolean discardFlag = false;

	private IntervalTreeMap<Boolean> intervalTreeMap=null;
	private IndexedBedReader bedReader =null;
	
	public VCFBedSetFilter()
		{
		}
		
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator r, VariantContextWriter w) {
		try {
			final VCFHeader h2=new VCFHeader(r.getHeader());
			addMetaData(h2);
			final VCFFilterHeaderLine filter = new VCFFilterHeaderLine(
					this.filterName,
					"Filtered with "+getProgramName()+", "+
					(this.inverse?" NOT  ":"")+
					"overlapping "+
					(this.tabixFile==null?this.treeMapFile:this.tabixFile)
					);
			
			
			if(!this.discardFlag) {
				h2.addMetaDataLine(filter);
			}
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(h2);
			w.writeHeader(h2);
			while(r.hasNext())
				{
				final VariantContext ctx= progress.watch(r.next());
				boolean set_filter=true;
				
				
				if(this.intervalTreeMap!=null) {
					if( this.intervalTreeMap.containsOverlapping(new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd())))
						{
						set_filter = false;	
						}
					}
				
				else 
					{
					final CloseableIterator<BedLine> iter = this.bedReader.iterator(
							ctx.getContig(),
							ctx.getStart()-1,
							ctx.getEnd()+1
							);
					while(iter.hasNext())
						{
						final BedLine bed = iter.next();
						if(!ctx.getContig().equals(bed.getContig())) continue;
						if(ctx.getStart() > bed.getEnd() ) continue;
						if(ctx.getEnd() < bed.getStart() ) continue;
						set_filter=false;
						break;
						}
					CloserUtil.close(iter);
					}
				
				if(this.inverse) set_filter=!set_filter;
				
				
				
				if(!set_filter)
					{
					w.add(ctx);
					continue;
					}
				
				if(!this.discardFlag)
					{
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.filter(filter.getID());
					w.add(vcb.make());
					}
				
				if(w.checkError()) break;
				}
			progress.finish();
			return RETURN_OK;
		} catch(Exception err) {
			LOG.error(err);
			return -1;
			}
		}

	
	@Override
	public int doWork(List<String> args) {
		try
			{
			if(this.tabixFile==null && this.treeMapFile==null)
				{
				LOG.error("Undefined tabix or memory file");
				return -1;
				}
			else if(this.tabixFile!=null && this.treeMapFile!=null)
				{
				LOG.error("You cannot use both options: treemap/tabix");
				return -1;
				}
			else if( this.tabixFile!=null) {
				LOG.info("opening Bed "+this.tabixFile);
				this.bedReader= new IndexedBedReader(this.tabixFile);
				}
			else 
				{
				LOG.info("opening Bed "+this.treeMapFile);
				this.intervalTreeMap  = super.readBedFileAsBooleanIntervalTreeMap(this.treeMapFile);
				}
			
			if(this.filterName==null || this.filterName.trim().isEmpty())
				{
				LOG.error("Undefined filter name");
				return -1;
				}
			
			return doVcfToVcf(args, outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.bedReader);
			this.bedReader = null;
			this.intervalTreeMap=null;
			}
		}
	
	public static void main(final String[] args) throws Exception
		{
		new VCFBedSetFilter().instanceMainWithExit(args);
		}
}
