/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
BEGIN_DOC
 


END_DOC
 */
@Program(name="vcfgnomad",
	description="Peek annotations from gnomad",
	keywords={"vcf","annotation","gnomad"},
	modificationDate="20200702"
)
public class VcfGnomad extends OnePassVcfLauncher {
	
	private static final Logger LOG = Logger.build(VcfGnomad.class).make();
	
	
	@Parameter(names={"-g","--gnomad"},description="Path to Gnomad VCF file.",required=true)
	private Path gnomadPath =null;
	@Parameter(names={"--bufferSize"},description="When we're looking for variant in Gnomad, load the variants for 'N' bases instead of doing a random access for each variant. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int gnomadBufferSize= 10_000;
	@Parameter(names={"--fields"},description="AF fields to peek-up from gnomad. Space/comma/semicolon separated")
	private String infoFieldStr="popmax,nfe";

	
	

	@Parameter(names={"--noUpdateId"},description="do Not Update ID if it is missing in user's variant")
	private boolean doNotUpdateId=false;
	@Parameter(names={"-gf","--gnomadFilter"},description="if defined, add this FILTER when any variant [CHROM:POS:REF] is found in nomad")
	private String inGnomadFilterName=null;
	@Parameter(names={"-of","--overlapFilter"},description="if defined, add this FILTER when any variant overlapping [CHROM:POS] is found in nomad")
	private String overlapGnomadFilterName=null;
	@Parameter(names={"-gnomad-filter-prefix","--gnomad-filter-prefix"},description="[20181214] if not empty, include the Gnomad FILTERs using this prefix.")
	private String filteredInGnomadFilterPrefix="GNOMAD";
	@Parameter(names={"--genome"},description="[20180327] For @MKarakachoff : genome only, don't use 'exome data'")
	private boolean useGenomeOnly = false;
	@Parameter(names={"--exclude"},description="[20180327] exclude gnomad INFO field matching this regular expression. Empty: accept all")
	private String excludePatternStr = "controls|non_cancer|non_neuro|non_topmed";
	@Parameter(names={"--ani"},description="[20190311] for allele numbers 'AN' to be variant-count-type=Integer (not 'A' as declared in gnomad). Deprecated : always true",hidden=true)
	private boolean alleleNumber_is_integer = true;
	@Parameter(names={"--ignore-error0"},description="[20190429] ignore error when gnomad/INFO is found twice for the same position. I found the error after a liftover to hg38. see https://twitter.com/yokofakun/status/1122814203381858305")
	private boolean ignore_info_found_twice = false;
	
	private BufferedVCFReader gnomadReader = null;
	private ContigNameConverter ctgNameConverter = null;
	private Set<String> gnomad_info_af_attributes = null;
	
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
		
		this.gnomad_info_af_attributes = Arrays.stream(this.infoFieldStr.split("[ ,;\t]+")).
				filter(S->!StringUtils.isBlank(S)).
				map(S->(S.startsWith("AF_")?S:"AF_"+S)).
				collect(Collectors.toSet());
		if(this.gnomad_info_af_attributes.isEmpty()) {
			LOG.error("No INFO attribute defined");
			return -1;
		}
		
		final VCFHeader gnomadHeader = this.gnomadReader.getHeader();
		for(final String field: gnomad_info_af_attributes) {
			final VCFInfoHeaderLine info = gnomadHeader.getInfoHeaderLine(field);
			if(info==null) {
				LOG.error("field INFO/"+field +" is undefined in "+this.gnomadPath);
				return -1;
				}
			if(info.getCountType()!=VCFHeaderLineCount.A) {
				LOG.warn("field INFO/"+field +" count-type is not 'A' but "+info.getCountType());
				}
			}
		
		/* do not keep those INFO in memory */
		final List<String> removeAtt =  this.gnomadReader.getHeader().getInfoHeaderLines().stream().
				map(H->H.getID()).
				filter(ID->!this.gnomad_info_af_attributes.contains(ID)).
				collect(Collectors.toList());
		
		this.gnomadReader.setSimplifier(V->{
			final VariantContextBuilder vcb = new VariantContextBuilder(V);
			vcb.rmAttributes(removeAtt);
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
		
		final VCFHeader gnomadHeader = this.gnomadReader.getHeader();
			
		if(!StringUtil.isBlank(this.filteredInGnomadFilterPrefix))
			{
			for(final VCFFilterHeaderLine fh: gnomadHeader.getFilterLines())
				{
				if(fh.getID().equals(VCFConstants.PASSES_FILTERS_v4)) continue;
				final String fid = this.filteredInGnomadFilterPrefix+"_"+ome.toUpperCase()+"_"+ fh.getID();
				final VCFFilterHeaderLine fh2 = new VCFFilterHeaderLine(
						fid,
						"[gnomad-"+ome+"]" + fh.getDescription()
						);
				h2.addMetaDataLine(fh2);
				}
			}
			
			
		
		
		if(!StringUtil.isBlank(this.inGnomadFilterName)) {
			h2.addMetaDataLine(new VCFFilterHeaderLine(this.inGnomadFilterName,"Variant CHROM/POS/REF was found in gnomad"));
			}
		
		if(!StringUtil.isBlank(this.overlapGnomadFilterName)) {
			h2.addMetaDataLine(new VCFFilterHeaderLine(this.overlapGnomadFilterName,"Gnomad Variant was found overlapping the variant, not necessarily at the same CHROM/POS"));
			}
		
		JVarkitVersion.getInstance().addMetaData(this, h2);
		out.writeHeader(h2);
		
		while(iter.hasNext()) {
			final VariantContext ctx = iter.next();
			
			final Set<String> filters = new HashSet<>(ctx.getFilters());

			
			
		
			if(ctx.getContig().equals("MT") || ctx.getContig().equals("chrM") || ctx.getContig().contains("_")) {
				out.add(ctx);
				continue;
				}
			
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			final List<Allele> alternateAlleles = ctx.getAlternateAlleles();
			String newid = null;
			boolean set_filter_ctx_is_in_gnomad = false;
			boolean found_gnomad_overlapping_variant = false;
			
			
			// variant overlapping 'ctx'
			final List<VariantContext> overlappingVariants = this.findOverlapping(ctx);
			if(!overlappingVariants.isEmpty()) found_gnomad_overlapping_variant = true;
			
			
			final List<VariantContext> gnomadVariants = overlappingVariants.
						stream().
						filter(V->V.getStart()==ctx.getStart() && V.getReference().equals(ctx.getReference())).
						collect(Collectors.toList());

				
				
				if(!gnomadVariants.isEmpty()) {
					set_filter_ctx_is_in_gnomad=true;
					
					// set new id ?
					if( newid == null) {
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
							map(F->this.filteredInGnomadFilterPrefix+"_"+ome.toUpperCase()+"_"+F).
							collect(Collectors.toList())
							);
						}
					}
				
				
				// loop over each field
				for(final String infoField: this.gnomad_info_af_attributes)
					{
					if(gnomadVariants.stream().noneMatch(V->V.hasAttribute(infoField))) continue;
					
					final double numbers[]=new double[alternateAlleles.size()];
					Arrays.fill(numbers,0.0);
					for(int x=0;x< alternateAlleles.size();++x)
						{
						final Allele alt=alternateAlleles.get(x);
						if(alt.equals(Allele.SPAN_DEL)) continue;
						for(final VariantContext gv:gnomadVariants)
							{
							int idx = gv.getAlternateAlleles().indexOf(alt);
							if(idx==-1) continue;
							if(!gv.hasAttribute(infoField)) continue;
							final List<Double> array = gv.getAttributeAsDoubleList(infoField,0.0);
							if(idx>=array.size()) continue;
							numbers[x] = array.get(idx);
							}
						}
					vcb.attribute(infoField, numbers);
					}
				
			if(set_filter_ctx_is_in_gnomad && !StringUtil.isBlank(this.inGnomadFilterName)) {
				filters.add(this.inGnomadFilterName);
				}
			
			if(found_gnomad_overlapping_variant && !StringUtil.isBlank(this.overlapGnomadFilterName)) {
				filters.add(this.overlapGnomadFilterName);
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
