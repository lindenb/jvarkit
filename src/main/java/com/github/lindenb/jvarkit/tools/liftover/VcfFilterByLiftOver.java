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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.liftover;

import java.io.File;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Function;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;
/**

```bash
```

END_DOC
*/
@Program(
		name="vcffilterbyliftover",
		description="Add FILTER(s) to a variant when it is known to map elsewhere after liftover.",
		keywords={"vcf","liftover"},
		modificationDate="20210603",
		creationDate="20190418"
		)
public class VcfFilterByLiftOver extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(VcfFilterByLiftOver.class).make();


	@Parameter(names={"-f","--chain"},description="LiftOver file.",required=true)
	private File liftOverFile = null;
	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = 1.0 ;
	@Parameter(names={"-d"},description="initial distance. See option -D.",converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_distance = 1_000;
	@Parameter(names={"-D"},description="final distance. We look for weird distance between the current and the previous variant on the same contig."
			+ "Two variants initially distance < d should have a distance <D after lift over. ",converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int max_distance = 1_500;
	@Parameter(names={"--no-validation"},description="Disable dictionary validation")
	private boolean disableValidation = false;
	
	@Parameter(names={"--lifted-filtered-vcf"},
		description="Another VCF in the destination reference. "
				+ "Input VCF will be filtered if this vcf is FILTERed at the same lifted position."
		)
	private Path anotherVcfFiltered = null;


	private VCFReader anotherVcfReader = null;
	
	
	private int distance(final Locatable L1,final Locatable L2) {
		if( CoordMath.overlaps(
				L1.getStart(), L1.getEnd(), 
				L2.getStart(), L2.getEnd())) return 0;
		return L1.getEnd() < L2.getStart()?
				L2.getStart()-L1.getEnd() : 
				L1.getStart()-L2.getEnd()
				;
		}
	
	private static String normalizeContig(final Locatable loc) {
		String str = loc.getContig().trim().toLowerCase();
		if(str.startsWith("chr")) str=str.substring(3);
		while(str.startsWith("0")) str=str.substring(1);
		if(str.equals("mt"))  str="m";
		return str;
		}
	
	private boolean sameContig(final Locatable a,final Locatable b) {
		return normalizeContig(a).equals(normalizeContig(b));
	}
	
	@Override
	protected int doVcfToVcf(final String inputName,final  VCFIterator in, final  VariantContextWriter out) {

		final ContigNameConverter anotherVcfCtgConverter;
		if(this.anotherVcfReader!=null) {
			final SAMSequenceDictionary dict =this.anotherVcfReader.getHeader().getSequenceDictionary();
			if(dict!=null) {
				anotherVcfCtgConverter = ContigNameConverter.fromOneDictionary(dict);
				}
			else
				{
				anotherVcfCtgConverter = ContigNameConverter.getIdentity();
				}
			}
		else
			{
			anotherVcfCtgConverter = ContigNameConverter.getIdentity();
			}
		
		final  LiftOver liftOver=new LiftOver(this.liftOverFile);
		liftOver.setLiftOverMinMatch(this.userMinMatch);

		
		final VCFHeader header = in.getHeader();
		final SAMSequenceDictionary dict = header.getSequenceDictionary();
		
		if(!this.disableValidation && dict!=null && !dict.isEmpty())  {
			liftOver.validateToSequences(dict);
		}
		
		final VCFHeader header2 = new VCFHeader(header);
		
		
		final VCFFilterHeaderLine filterLiftOverFailed = new VCFFilterHeaderLine(
				"LIFTOVER_FAILED", "liftover failed "+this.liftOverFile);
		header2.addMetaDataLine(filterLiftOverFailed);

		
		final VCFFilterHeaderLine filterNoSameContig = new VCFFilterHeaderLine(
				"LIFTOVER_OTHER_CTG", "Variant is mapped to another contig after liftover with "+this.liftOverFile);
		header2.addMetaDataLine(filterNoSameContig);

		final VCFInfoHeaderLine infoLiftOverPos = new VCFInfoHeaderLine(
				"LIFTOVER_LOC", 
				1,
				VCFHeaderLineType.String,
				"Position of the variant liftover-ed with "+this.liftOverFile
				);
		header2.addMetaDataLine(infoLiftOverPos);

		final VCFFilterHeaderLine filterDistantFromPrev = new VCFFilterHeaderLine(
				"LIFTOVER_DISTANT", "After liftover the distance (< "+min_distance+") with the previous variant is unusual( > "+max_distance+") after liftover with "+this.liftOverFile);
		header2.addMetaDataLine(filterDistantFromPrev);
		
		/* transfert filters to new header */
		if(this.anotherVcfReader!=null)  {
			this.anotherVcfReader.
				getHeader().
				getFilterLines().
				forEach(F->header2.addMetaDataLine(F));
			}
		
		
		JVarkitVersion.getInstance().addMetaData(this, header2);
		out.writeHeader(header2);
		
		Locatable prevCtx=null;
		Locatable prevLifted=null;
		
		
		final Function<Locatable, String> interval2str = R->R.getContig()+"|"+R.getStart()+"|"+R.getEnd();
		
		while(in.hasNext())
			{
			final VariantContext ctx= in.next();
			
			if(prevCtx!=null && !prevCtx.getContig().equals(ctx.getContig())) {
				prevCtx = null;
				prevLifted = null;
				}
			
			final Interval lifted=liftOver.liftOver(
					new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd(),
					false,//negative strand
					interval2str.apply(ctx)
					));
			// lifover failed
			if(lifted==null) {
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				vcb.filter(filterLiftOverFailed.getID());
				out.add(vcb.make());
				}
			// another contig
			else if(lifted!=null && !sameContig(lifted,ctx))
				{
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				
				vcb.filter(filterNoSameContig.getID());
				vcb.attribute(infoLiftOverPos.getID(),interval2str.apply(lifted));

				
				out.add(vcb.make());
				}
			// strange distance
			else if(prevCtx!=null && lifted!=null && prevLifted!=null &&
					prevCtx.getContig().equals(ctx.getContig()) &&
					sameContig(prevLifted,lifted) &&
					distance(prevCtx, ctx) < this.min_distance &&
					distance(prevLifted, lifted) > this.max_distance
					)
				{
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				vcb.filter(filterDistantFromPrev.getID());
				vcb.attribute(infoLiftOverPos.getID(),interval2str.apply(lifted));
				out.add(vcb.make());
				}
			//filtered in anotherVcf
			else if(this.anotherVcfReader!=null) {
				final Set<String> found_filtered = new HashSet<>();
				final String ctg2 = anotherVcfCtgConverter.apply(lifted.getContig());
				if(!StringUtils.isBlank(ctg2)) {
					try(CloseableIterator<VariantContext> iter2= this.anotherVcfReader.query(ctg2, lifted.getStart(),lifted.getEnd())) {
						while(iter2.hasNext()) {
							final VariantContext ctx2 = iter2.next();
							if(!ctx2.isFiltered() || ctx2.getLengthOnReference()!=ctx.getLengthOnReference()) continue;
							found_filtered.addAll(ctx2.getFilters());
							break;
							}
						}
					}
				if(found_filtered.isEmpty()) {
					out.add(ctx);
					}
				else
					{
					//add previous
					found_filtered.addAll(ctx.getFilters());
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					vcb.filters(found_filtered);
					out.add(vcb.make());
					}
				}
			else
				{
				out.add(ctx);
				}
			prevCtx = ctx;
			prevLifted = lifted;
			}
		return 0;
		}

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeVcf() {
		if(this.liftOverFile==null)
			{
			LOG.error("LiftOver file is undefined.");
			return -1;
			}
		IOUtil.assertFileIsReadable(this.liftOverFile);
		if(this.anotherVcfFiltered!=null) {
			this.anotherVcfReader = VCFReaderFactory.
					makeDefault().
					open(this.anotherVcfFiltered,true);
			}
		return super.beforeVcf();
		}
	@Override
	protected void afterVcf() {
		CloserUtil.close(this.anotherVcfReader);
		super.afterVcf();
		}

	public static void main(final String[] args) {
		new VcfFilterByLiftOver().instanceMainWithExit(args);
		}
	}
