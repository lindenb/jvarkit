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
package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;


/**

BEGIN_DOC

### Examples

```
$java -jar dist/vcfbedsetfilter.jar -f MYFILTER - -B in.bed in.vcf 
```

## history:

2191104: changed the logic which was wrongly defined in the documentation. :-/

END_DOC
*/



@Program(name="vcfbedsetfilter",
	description="Set FILTER for VCF if intersects with BED.",
	keywords={"vcf","bed","filter"},
	modificationDate="20210417",
	creationDate="20150415",
	biostars={9465226}
	)
public class VCFBedSetFilter extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VCFBedSetFilter.class).make();


	@Parameter(names={"-f","--filter"},description="FILTER name. If `--filter` is empty, FILTERED variant will be discarded.")
	private String filterName = "VCFBED";
	@Parameter(names={"-e","--exclude","--blacklist"},description="Tribble or Tabix bed file containing the regions to be FILTERED. Must be indexed with tribble or tabix, or use '--fast' to load in memory.")
	private File tabixBlackFile = null;
	@Parameter(names={"-i","--include","--whitelist"},description="Tribble or Tabix bed file containing the regions to be accepted. Regions NOT overlapping those regions will be FILTERED. Must be indexed with tribble or tabix, or use '--fast' to load in memory.")
	private File tabixWhiteFile = null;
	@Parameter(names={"--fast","--memory"},description="Load the bed in memory: faster than tribble/tabix but memory consumming)")
	private boolean useInMemory = false;
	@Parameter(names={"-x","--extend"},description="Extend the variant coordinates per 'x' bases. " + DistanceParser.OPT_DESCRIPTION ,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int extend_bases = 0;
	@Parameter(names={"--debug"},description="debug what's happening",hidden=true)
	private boolean debug = false;
	@Parameter(names={"--min-bed-fraction"},description="Min BED fraction overlap after extension. Only consider BED records if variant overlap >= 'x' percent of bed length. "+FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private Double min_bed_fraction_overlap = null;
	@Parameter(names={"--min-vc-fraction"},description="Min Variant fraction overlap after extension. Only consider BED records if bed overlap >= 'x' percent of vc length. "+ FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private Double min_vc_fraction_overlap = null;

	

	private IntervalTreeMap<Interval> intervalTreeMap=null;
	private IndexedBedReader bedReader =null;
	private File tabixFile = null;
	
	public VCFBedSetFilter()
		{
		}
	
	private boolean testOverlap(final int ctx_start,final int ctx_end,final Locatable bed) {
		if(!CoordMath.overlaps(bed.getStart(), bed.getEnd(),ctx_start,ctx_end)) return false;
		// get overlap bounds
		final int x1 = Math.max(bed.getStart(), ctx_start); 
		final int x2 = Math.min(bed.getEnd(), ctx_end); 
		final double len_x =  CoordMath.getLength(x1, x2);

		
		if(min_bed_fraction_overlap!=null) {
			final double len_bed = bed.getLengthOnReference();
			if(len_x/len_bed < this.min_bed_fraction_overlap) return false;
		}
		if(min_vc_fraction_overlap!=null) {
			final double len_vc =  CoordMath.getLength(ctx_start, ctx_end);
			if(len_x/len_vc < this.min_vc_fraction_overlap) return false;
		}
		return true;
	}
	
	
	@Override
	protected int beforeVcf() {
		if(this.tabixBlackFile==null && this.tabixWhiteFile==null) {
			LOG.error("include / exclude file are both undefined");
			return -1;
		}
		if(this.tabixBlackFile!=null && this.tabixWhiteFile!=null) {
			LOG.error("include / exclude file are both defined");
			return -1;
		}
		
		this.tabixFile = (this.tabixBlackFile==null?this.tabixWhiteFile:this.tabixBlackFile);
		
		try {
			if(this.useInMemory) {
				this.intervalTreeMap  = new IntervalTreeMap<>();
				final BedLineCodec bedCodec=new BedLineCodec();
	
				try(BufferedReader br= IOUtils.openFileForBufferedReading(this.tabixFile)) {
					br.lines().
						filter(line->!StringUtil.isBlank(line)).
						filter(line->!BedLine.isBedHeader(line)).
						map(line->bedCodec.decode(line)).
						filter(B->B!=null).
						map(B->B.toInterval()).
						forEach(L->intervalTreeMap.put(L,L));
					}
				}
			else 
				{
				this.bedReader= new IndexedBedReader(this.tabixFile);
				}
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		
		return 0;
		}
	
	@Override
	protected void afterVcf() {
		this.intervalTreeMap=null;
		CloserUtil.close(this.bedReader);
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator r, final VariantContextWriter w) {
			final Set<String> contigs_not_found = new HashSet<>();
			final VCFHeader h2=new VCFHeader(r.getHeader());
						
			final VCFFilterHeaderLine filter;
			
			if(!StringUtil.isBlank(this.filterName)) {
				filter = new VCFFilterHeaderLine(
					this.filterName,
					"Variant "+
					(this.tabixBlackFile!=null?" overlapping any bed record of ":" that do not overlap any bed record of ")+
					tabixFile
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
			final SAMSequenceDictionary dict = r.getHeader().getSequenceDictionary();
			w.writeHeader(h2);
			while(r.hasNext())
				{
				final VariantContext ctx = r.next();
				boolean set_filter = false;
				final String convert_contig = ctgNameConverter.apply(ctx.getContig());
				final int ctx_start = Math.max(1, ctx.getStart() - this.extend_bases);
				final int ctx_end;
				if(dict==null) {
					ctx_end = ctx.getEnd() + this.extend_bases ;
					}
				else
					{
					final SAMSequenceRecord ssr = dict.getSequence(ctx.getContig());
					if(ssr==null) throw new JvarkitException.ContigNotFoundInDictionary(ctx.getContig(), dict);
					ctx_end = Math.min(ssr.getSequenceLength(),ctx.getEnd() + this.extend_bases);
					}
				
				if(StringUtil.isBlank(convert_contig))
					{
					if(debug) LOG.warn("Cannot convert contig "+ ctx.getContig() );
					if(contigs_not_found.size()<100) {
						if(contigs_not_found.add(ctx.getContig()))
							{
							LOG.warn("Cannot convert variant contig "+ctx.getContig()+" to bed file. (Contig is not in BED file)");
							}
						}
					set_filter = false;
					}
				else if(this.intervalTreeMap!=null) {
					if( this.intervalTreeMap.getOverlapping(new SimpleInterval(convert_contig,ctx_start,ctx_end)).
							stream().
							anyMatch(LOC->testOverlap(ctx_start, ctx_end, LOC)))
						{
						if(debug) LOG.warn("treemap.overlap=true set Filter=true for "+ctx.getContig()+":"+ctx.getStart() );
						set_filter = true;	
						}
					}
				
				else 
					{
					try(CloseableIterator<BedLine> iter = this.bedReader.iterator(
							convert_contig,
							ctx_start-1,
							ctx_end+1
							)) {
						while(iter.hasNext())
							{
							final BedLine bed = iter.next();
							if(bed==null ) continue;
							if(!testOverlap(ctx_start, ctx_end, bed)) continue;
							if(debug) LOG.warn("tabix=true set Filter=true for "+ctx.getContig()+":"+ctx.getStart() );
							set_filter=true;
							break;
							}
						} catch(IOException err) {
							throw new RuntimeIOException(err);
						}
					}
				
				/* variant is in whitelist, we want to keep it */
				if(this.tabixWhiteFile!=null) {
					set_filter=!set_filter;
					if(debug) LOG.warn("inverse. Now Filter="+set_filter+" for "+ctx.getContig()+":"+ctx.getStart() );
					}
				
				
				if(!set_filter)
					{
					if(debug) LOG.warn("no filter. Writing "+ctx.getContig()+":"+ctx.getStart() );
					w.add(ctx);
					continue;
					}
				
				if(filter != null )
					{
					if(debug) LOG.warn("adding filter. Writing "+ctx.getContig()+":"+ctx.getStart() );
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.filter(filter.getID());
					w.add(vcb.make());
					}
				else
					{
					if(debug) LOG.warn("Ignoring "+ctx.getContig()+":"+ctx.getStart() );
					}
				}
			if(!contigs_not_found.isEmpty()) {
				LOG.warn(
					"The following contigs were not found: "+
					String.join(" ", contigs_not_found)+ 
					"..."
					);
				}
			return 0;
			}
	
	public static void main(final String[] args) throws Exception
		{
		new VCFBedSetFilter().instanceMainWithExit(args);
		}
	}
