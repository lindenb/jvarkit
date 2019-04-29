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
package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;


/**

BEGIN_DOC

### Examples

```
$java -jar dist/vcfbedsetfilter.jar -f MYFILTER - -B in.bed in.vcf 
```

END_DOC
*/



@Program(name="vcfbedsetfilter",
	description="Set FILTER for VCF if intersects with BED.",
	keywords={"vcf","bed","filter"},
	modificationDate="20190426"
	)
public class VCFBedSetFilter extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFBedSetFilter.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-f","--filter"},description="FILTER name. "
			+ "Filter is **set** if the variant overlaps any BED region, unless `--inverse` is set. "
			+ "If `--filter` is empty, FILTERED variant will be discarded.")
	private String filterName = "VCFBED";

	@Parameter(names={"-i","--inverse"},description="Inverse selection: FILTER will be **set** for a Variant overlaping no bed record. Variant overlapping any bed record remains unfiltered.")
	private boolean inverse = false;

	@Parameter(names={"-B","--bed","-m","--map"},description="Tribble or Tabix bed file. Must be indexed with tribble or tabix, or use '--fast' to load in memory.", required=true)
	private File tabixFile = null;

	@Parameter(names={"--fast","--memory"},description="Load the bed in memory: faster than tribble/tabix but memory consumming)")
	private boolean useInMemory = false;

	@Parameter(names={"-x","--extend"},description="Extend the variant coordinates per 'x' bases. " + DistanceParser.OPT_DESCRIPTION ,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int extend_bases = 0;


	private IntervalTreeMap<Boolean> intervalTreeMap=null;
	private IndexedBedReader bedReader =null;
	
	public VCFBedSetFilter()
		{
		}
		
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator r,
			final VariantContextWriter w) {
		try {
			final Set<String> contigs_not_found = new HashSet<>();
			final VCFHeader h2=new VCFHeader(r.getHeader());
			
			
			
			final VCFFilterHeaderLine filter;
			
			if(!StringUtil.isBlank(this.filterName)) {
				filter = new VCFFilterHeaderLine(
			
					this.filterName,
					"Variant "+
					(this.inverse?" NOT  ":"")+
					"overlapping any bed record of "+
					this.tabixFile
					);
				h2.addMetaDataLine(filter);
				}
			else
				{
				filter = null;
				}
			JVarkitVersion.getInstance().addMetaData(this, h2);
			
			final ContigNameConverter ctgNameConverter;
			if(this.bedReader!=null)
				{
				ctgNameConverter = ContigNameConverter.fromContigSet(this.bedReader.getContigs());
				}
			else
				{
				ctgNameConverter  = ContigNameConverter.fromIntervalTreeMap(this.intervalTreeMap);
				}
			
			
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(h2).logger(LOG).build();
			w.writeHeader(h2);
			while(r.hasNext())
				{
				final VariantContext ctx= progress.apply(r.next());
				boolean set_filter = true;
				final String convert_contig = ctgNameConverter.apply(ctx.getContig());
				final int ctx_start = Math.max(1, ctx.getStart() - this.extend_bases);
				final int ctx_end =  ctx.getEnd() + this.extend_bases ;
				
				if(StringUtil.isBlank(convert_contig))
					{
					if(contigs_not_found.size()<100) {
						if(contigs_not_found.add(ctx.getContig()))
							{
							LOG.warn("Cannot convert variant contig "+ctx.getContig()+" to bed file. (Contig is not in BED file)");
							}
						}
					set_filter = false;
					}
				else if(this.intervalTreeMap!=null) {
					if( this.intervalTreeMap.containsOverlapping(new Interval(convert_contig,ctx_start,ctx_end)))
						{
						set_filter = false;	
						}
					}
				
				else 
					{
					try(CloseableIterator<BedLine> iter = this.bedReader.iterator(
							convert_contig,
							ctx.getStart()-1,
							ctx.getEnd()+1
							)) {
						while(iter.hasNext())
							{
							final BedLine bed = iter.next();
							if(bed==null || !CoordMath.overlaps(bed.getStart(), bed.getEnd(),ctx_start,ctx_end)) continue;
							set_filter=false;
							break;
							}
						}
					}
				
				
				if(this.inverse) {
					set_filter=!set_filter;
				}
				
				
				if(!set_filter)
					{
					w.add(ctx);
					continue;
					}
				
				if(filter != null )
					{
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.filter(filter.getID());
					w.add(vcb.make());
					}
				}
			progress.close();
			if(!contigs_not_found.isEmpty()) {
				LOG.warn(
					"The following contigs were not found: "+
					String.join(" ", contigs_not_found)+ 
					"..."
					);
				}
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			if(this.tabixFile==null)
				{
				LOG.error("Undefined tabix or memory file");
				return -1;
				}
			
			if(this.useInMemory) {
				this.intervalTreeMap  = super.readBedFileAsBooleanIntervalTreeMap(this.tabixFile);
				}
			else 
				{
				this.bedReader= new IndexedBedReader(this.tabixFile);
				}
			
			return doVcfToVcf(args, this.outputFile);
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
