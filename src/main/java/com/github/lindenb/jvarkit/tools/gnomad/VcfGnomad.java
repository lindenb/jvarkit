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
package com.github.lindenb.jvarkit.tools.gnomad;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC
 
The tool was redesigned on July 2nd, 2020.

## Example:

```
java -jar dist/vcfgnomad.jar -g src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz src/test/resources/test_vcf01.vcf
```

END_DOC
 */
@Program(name="vcfgnomad",
	description="Peek annotations from gnomad",
	keywords={"vcf","annotation","gnomad"},
	modificationDate="20200702",
	creationDate="20170407"
)
public class VcfGnomad extends OnePassVcfLauncher {
	
	private static final Logger LOG = Logger.build(VcfGnomad.class).make();
	
	
	@Parameter(names={"-g","--gnomad"},description="Path to Indexed Gnomad VCF file.",required=true)
	private Path gnomadPath =null;
	@Parameter(names={"--bufferSize"},description= BufferedVCFReader.OPT_BUFFER_DESC+" "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int gnomadBufferSize= 10_000;
	@Parameter(names={"-F","--fields"},description="AF fields to peek-up from gnomad. Space/comma/semicolon separated")
	private String infoFieldStr="AF_popmax,AF_nfe";
	@Parameter(names={"--noUpdateId"},description="do Not Update ID if it is missing in user's variant")
	private boolean doNotUpdateId=false;
	@Parameter(names={"--prefix"},description="If not empty, include the Gnomad FILTERs using this prefix.")
	private String filteredInGnomadFilterPrefix="GNOMAD";
	
	@Parameter(names={"--min-af"},description="Min allele frequency",converter=FractionConverter.class,splitter=NoSplitter.class)
	private double min_af = 0.0;
	@Parameter(names={"--max-af"},description="Max allele frequency",converter=FractionConverter.class,splitter=NoSplitter.class)
	private double max_af = 1.0;
	@Parameter(names={"--debug"},description="debug",hidden=true)
	private boolean debug = false;

	
	private BufferedVCFReader gnomadReader = null;
	private ContigNameConverter ctgNameConverter = null;
	private final Set<String> gnomad_info_af_attributes = new HashSet<>();
	
	private String toString(final VariantContext vc) {
	return vc.getContig()+":"+vc.getStart()+":"+vc.getReference();
	}

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeVcf() {
		try {
			final VCFReader r = VCFReaderFactory.makeDefault().open(this.gnomadPath,true);
			this.gnomadReader = new BufferedVCFReader(r, this.gnomadBufferSize);
			this.ctgNameConverter = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(r.getHeader()));
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		
		
		final VCFHeader gnomadHeader = this.gnomadReader.getHeader();
		for(final String field: this.infoFieldStr.split("[ ,;\t]+")) {
			if(StringUtils.isBlank(field)) continue;
			
			final VCFInfoHeaderLine info = gnomadHeader.getInfoHeaderLines().
					stream().
					filter(F->F.getID().equalsIgnoreCase(field)).
					findFirst().orElse(null);
			if(info==null) {
				LOG.error("field INFO/"+field +" is undefined in "+this.gnomadPath);
				return -1;
				}
			if(!field.equals(info.getID())) {
				LOG.warn("changed user field INFO/"+field +" to INFO/"+info.getID());
				}
			if(info.getCountType()!=VCFHeaderLineCount.A) {
				LOG.warn("field INFO/"+field +" count-type is not 'A' but "+info.getCountType());
				}
			this.gnomad_info_af_attributes.add(info.getID());
			}
		if(this.gnomad_info_af_attributes.isEmpty()) {
			LOG.error("No INFO attribute defined");
			return -1;
			}

		
		/* do not keep those INFO in memory */
		final List<String> removeAtt =  this.gnomadReader.getHeader().getInfoHeaderLines().stream().
				map(H->H.getID()).
				filter(ID->!this.gnomad_info_af_attributes.contains(ID)).
				collect(Collectors.toList());
		
		this.gnomadReader.setSimplifier(V->{
			final VariantContextBuilder vcb = new VariantContextBuilder(V);
			vcb.rmAttributes(removeAtt).noGenotypes();
			return vcb.make();
			});
		
		return 0;
		}
	
	@Override
	protected void afterVcf() {
		try {
			this.gnomadReader.close();
			}
		catch(final Throwable err) {
			LOG.error(err);
			}
		}
	
	
	/** find matching variant in tabix file, use a buffer to avoid multiple random accesses */
	private List<VariantContext> findOverlapping(final VariantContext userVariantCtx)
		{
		final String normContig = this.ctgNameConverter.apply(userVariantCtx.getContig());
		final Locatable loc;
		if(StringUtil.isBlank(normContig)) return Collections.emptyList();
		if(normContig.equals(userVariantCtx.getContig())) {
			loc = userVariantCtx;
			}
		else
			{
			loc = new SimpleInterval(normContig, userVariantCtx.getStart(), userVariantCtx.getEnd());
			}
		
		final List<VariantContext> list = new ArrayList<>();
		
		try( CloseableIterator<VariantContext> iter= this.gnomadReader.query(loc)) {
			while(iter.hasNext())
				{
				final VariantContext ctx = iter.next();
				list.add(ctx);
				}
			}
		return list;
		}
	
	
	
	
		
		
		
		
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter out
			)
		{
		final VCFHeader h0 = iter.getHeader();
		final VCFHeader h2 = new VCFHeader(h0);
		final String ome;
		if(this.gnomadPath.getFileName().toString().contains("genomes")) {
			ome= "genome";
			}
		else if(this.gnomadPath.getFileName().toString().contains("exomes")) {
			ome= "exome";
			}
		else
			{
			LOG.info("Cannot identify if gnomad is exome or genome :"+this.gnomadPath);
			return -1;
			}
		/* output FILTER name */
		final UnaryOperator<String> toNewFilter = S-> this.filteredInGnomadFilterPrefix+"_"+ome.toUpperCase()+"_"+S;
		/* output INFO name */
		final UnaryOperator<String> toNewInfo = S-> "gnomad_" + ome + "_"+ S.toUpperCase();
		
		
		final VCFHeader gnomadHeader = this.gnomadReader.getHeader();
		final VCFFilterHeaderLine filterFrequencyHeader;
		
		/* peek FILTER from GNOMAD */
		if(!StringUtil.isBlank(this.filteredInGnomadFilterPrefix))
			{
			for(final VCFFilterHeaderLine fh: gnomadHeader.getFilterLines())
				{
				if(fh.getID().equals(VCFConstants.PASSES_FILTERS_v4)) continue;
				final VCFFilterHeaderLine fh2 = new VCFFilterHeaderLine(
						toNewFilter.apply(fh.getID()),
						"[gnomad-"+ome+"]" + fh.getDescription()+" in "+ gnomadPath
						);
				if(this.debug) LOG.debug("adding vcf filter "+fh2.getID());
				h2.addMetaDataLine(fh2);
				}
			filterFrequencyHeader = new VCFFilterHeaderLine(this.filteredInGnomadFilterPrefix+"_"+ome.toUpperCase()+"_BAD_AF",
					"AF if not between "+this.min_af+"<= 'af' <="+this.max_af+" for "+this.gnomadPath
					);
			h2.addMetaDataLine(filterFrequencyHeader);
			if(this.debug) LOG.debug("adding filter "+ filterFrequencyHeader.getID());

			}
		else
			{
			filterFrequencyHeader = null;
			}
		
		/* peek INFO from GNOMAD */
		for(final String field: this.gnomad_info_af_attributes)
			{
			final VCFInfoHeaderLine fh = Objects.requireNonNull(gnomadHeader.getInfoHeaderLine(field));
			
			
			final VCFInfoHeaderLine fh2 = new VCFInfoHeaderLine(
					toNewInfo.apply(fh.getID()),
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"[gnomad-"+ome+"]" + fh.getDescription()+" in "+ gnomadPath
					);
			h2.addMetaDataLine(fh2);
			if(this.debug) LOG.debug("Adding INFO vcf "+fh2.getID());
			}
		
		
		final VCFInfoHeaderLine infoFlagContigStartRef = new VCFInfoHeaderLine(
				"IN_GNOMAD_"+ome.toUpperCase(),
				1,
				VCFHeaderLineType.Flag,
				"Variant CHROM/POS/REF was found in gnomad "+this.gnomadPath
				);
		h2.addMetaDataLine(infoFlagContigStartRef);
		if(this.debug) LOG.debug("adding vcfheader "+infoFlagContigStartRef.getID());

		final VCFInfoHeaderLine infoNumOverlapping = new VCFInfoHeaderLine(
				"N_GNOMAD_"+ome.toUpperCase(),
				1,
				VCFHeaderLineType.Integer,
				"Count Gnomad Variants that were found overlapping the user variant, not necessarily at the same CHROM/POS in "+this.gnomadPath
				);
		h2.addMetaDataLine(infoNumOverlapping);
		if(this.debug) LOG.debug("adding vcfheader "+infoNumOverlapping.getID());
		
		
		JVarkitVersion.getInstance().addMetaData(this, h2);
		out.writeHeader(h2);
		
		while(iter.hasNext()) {
			final VariantContext ctx = iter.next();
			
			final Set<String> filters = new HashSet<>(ctx.getFilters());

			
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			
			for(final String infoField: this.gnomad_info_af_attributes) {
				vcb.rmAttribute(toNewInfo.apply(infoField));
				}
			
			final List<Allele> alternateAlleles = ctx.getAlternateAlleles();
			String newid = null;
			
			
			// variant overlapping 'ctx'
			final List<VariantContext> overlappingVariants = this.findOverlapping(ctx);
			vcb.attribute(infoNumOverlapping.getID(), overlappingVariants.size());			

			
			final List<VariantContext> gnomadVariants = overlappingVariants.
						stream().
						filter(V->V.getStart()==ctx.getStart() && V.getReference().equals(ctx.getReference())).
						collect(Collectors.toList());

			double ctx_min_AF = 1.0;
			if(!ctx.isVariant() ||
				gnomadVariants.isEmpty() ||
				(alternateAlleles.size()==1 && alternateAlleles.get(0).equals(Allele.SPAN_DEL)
				)) {
				ctx_min_AF = 0.0;// not in gnomad
				for(final String infoField: this.gnomad_info_af_attributes) {
					final double numbers[]=new double[alternateAlleles.size()];
					Arrays.fill(numbers,0.0);
					vcb.attribute(toNewInfo.apply(infoField), numbers);
					}
				}
			else
				{
				vcb.attribute(infoFlagContigStartRef.getID(), true);
				
				// set new id ?
				if( newid == null && !ctx.hasID()) {
					newid = gnomadVariants.
						stream().
						filter(V->V.hasID()).
						map(V->V.getID()).
						findFirst().
						orElse(null);
					}
					
				// add FILTER(s)
				if(!StringUtil.isBlank(this.filteredInGnomadFilterPrefix)) {
					filters.addAll(
						gnomadVariants.
						stream().
						filter(V->V.isFiltered()).
						flatMap(V->V.getFilters().stream()).
						filter(F->!F.equals(VCFConstants.PASSES_FILTERS_v4)).
						map(F->toNewFilter.apply(F)).
						collect(Collectors.toList())
						);
					}
				
				// loop over each field
				for(final String infoField: this.gnomad_info_af_attributes) {
					final double numbers[]=new double[alternateAlleles.size()];
					Arrays.fill(numbers,0.0);
					for(int x=0;x< alternateAlleles.size();++x)
						{
						final Allele alt= alternateAlleles.get(x);
						if(alt.equals(Allele.SPAN_DEL)) continue;
						Double alt_af = null;
						for(final VariantContext gv:gnomadVariants)
							{
							final int idx = gv.getAlternateAlleles().indexOf(alt);
							// this ALT is NOT in gnomad variant
							if(idx==-1) {
								if(this.debug) LOG.debug(toString(ctx)+ " no allele "+alt+ " in gnomad "+ toString(gv));
								continue;
								}
							if(!gv.hasAttribute(infoField)) {
								if(this.debug) LOG.debug(toString(ctx)+ " no info "+ infoField +" in gnomad "+ toString(gv));
								continue;
								}
							final List<Double> array = gv.getAttributeAsDoubleList(infoField,0.0);
							if(idx>=array.size()) {
								if(this.debug) LOG.debug(toString(ctx)+ " array out of bound for "+ infoField + " in gnomad "+ toString(gv));
								continue;
								}
							alt_af = array.get(idx);
							if(this.debug) LOG.debug(toString(ctx)+ " "+infoField+" got AF="+alt_af);
							break;
							}
						if(alt_af==null) alt_af = 0.0;
						ctx_min_AF =  Math.min(ctx_min_AF  , alt_af);
						if(this.debug) LOG.debug(toString(ctx)+ " "+infoField+":"+alt_af+" min"+ctx_min_AF);
						numbers[x] = alt_af;
						}
					vcb.attribute(toNewInfo.apply(infoField), numbers);
					}
				}
			// test for frequency
			if(!(this.min_af<=ctx_min_AF && ctx_min_AF<=this.max_af) ) {
				if(this.debug) LOG.debug(toString(ctx)+ "fails  "+ this.min_af +"<="+ctx_min_AF+"<="+this.max_af);
				// skip variants
				if(filterFrequencyHeader==null) {
					continue;
					}
				else
					{
					if(this.debug) LOG.debug(toString(ctx)+ "set filter "+filterFrequencyHeader.getID());
					filters.add(filterFrequencyHeader.getID());
					}
				}
			
			if(!this.doNotUpdateId && !ctx.hasID() && !StringUtil.isBlank(newid)) vcb.id(newid);
			vcb.filters(filters);
			out.add(vcb.make());
			}
		return 0;
		}
	
	
	

public static void main(final String[] args) {
	new VcfGnomad().instanceMainWithExit(args);
	}
}
