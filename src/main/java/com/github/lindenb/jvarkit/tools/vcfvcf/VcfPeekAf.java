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
package com.github.lindenb.jvarkit.tools.vcfvcf;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;
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
/**
BEGIN_DOC
ALLELE_FREQUENCY_KEY
END_DOC

 */
@Program(name="vcfpeekaf",
		description="Peek the AF from another VCF",
		keywords={"vcf","annotation","af"},
		creationDate="20200624",
		modificationDate="20200904"
		)
public class VcfPeekAf extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VcfPeekAf.class).make();

	@Parameter(names={"-F","--database","--tabix","--resource"},description="An indexed VCF file. Source of the annotations",required=true)
	private Path resourceVcfFile = null;
	@Parameter(names={"-T","--tag"},description="INFO tag to put found frequency. empty: no extra tag.")
	private String frequencyTag = "";
	@Parameter(names={"-f","--filter"},description="soft FILTER the variant of this data if AF is not found or it greater > max-af or lower than min-af. If empty, just DISCARD the variant")
	private String filterStr = "";
	@Parameter(names={"-b","--buffer-size"},converter=DistanceParser.StringConverter.class, description=BufferedVCFReader.OPT_BUFFER_DESC+" "+DistanceParser.OPT_DESCRIPTION,splitter=NoSplitter.class)
	private int buffer_size = 10_000;
	@Parameter(names={"-l","--list"},description="List available AF peekers and exit.",help=true)
	private boolean list_peekers = false;
	@Parameter(names={"-p","--peeker"},description="AF Peeker name. Use option --list to get a list of peekers.",required=true)
	private String peekerName = null;
	@Parameter(names={"-t","--treshold","--max-af"},description="AF max treshold. Variant is accepted is computed AF <= treshold.",converter=FractionConverter.class,required=true)
	private double af_maximum = 1.0;
	@Parameter(names={"--min-af"},description="AF min treshold. Variant is accepted is computed AF >= treshold.",converter=FractionConverter.class,required=false)
	private double af_minimum = 0.0;
	@Parameter(names={"--no-alt"},description="Do not look at the alternate alleles concordance")
	private boolean disable_alt_concordance = false;
	@Parameter(names={"-P","--peek-info"},description="Name of INFO tag in the vcf database to extract the AF value for exractor .'Custom'"  )
	private String custom_peek_info_name = null;
	@Parameter(names={"--peek-id"},description="Peek database variant ID if it is missing in the processed VCF."  )
	private boolean peek_variant_id=false;

	
	/* a class extracting the allele frequency from another VCF */
	private abstract class AFPeeker
		{
		abstract String getName();
		abstract String getDescription();
		abstract void initialize(final VCFHeader h);
		@Override
		public String toString() {
			return getName();
			}
		/** compare ctx with other variants found in database, ignoring ALT status */
		abstract VariantContext applyAlt(final VariantContext ctx,final List<VariantContext> overlappers);
		/** compare ctx with other variants found in database, considering ALT status */
		abstract VariantContext applyIgnoringAlt(final VariantContext ctx,final List<VariantContext> overlappers);
		
		final VariantContext apply(final VariantContext ctx,final List<VariantContext> overlappers) {
			return disable_alt_concordance ? applyIgnoringAlt(ctx,overlappers):applyAlt(ctx, overlappers);
			}

		/** simplify database variant */
		abstract VariantContext sanitize(final VariantContext ctx);
		
		/** change variant in VCF, applying filters if needed . Return null if the variant is discarded. */
		VariantContext addFilters3(final VariantContext ctx,final VariantContextBuilder vcb,boolean accept) {
			final Set<String> old_filters = new HashSet<>(ctx.getFilters());
			if(!StringUtils.isBlank(VcfPeekAf.this.filterStr)) old_filters.remove(VcfPeekAf.this.filterStr);

			
			if(accept) {
				if(old_filters.isEmpty()) {
					vcb.passFilters();
					}
				else
					{
					vcb.filters(old_filters);
					}
				return vcb.make();
				}
			else
				{
				if(StringUtils.isBlank(VcfPeekAf.this.filterStr)) return null;
				old_filters.add(VcfPeekAf.this.filterStr);
				vcb.filters(old_filters);
				return vcb.make();
				}
			}

		
		VariantContext addFiltersIgnoreAlt(final VariantContext ctx,final OptionalDouble optFreq) {
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			final boolean foundLt= !optFreq.isPresent() || (optFreq.getAsDouble() >= VcfPeekAf.this.af_minimum && optFreq.getAsDouble() <= VcfPeekAf.this.af_maximum);
			if(optFreq.isPresent() && !StringUtils.isBlank(VcfPeekAf.this.frequencyTag)) {
				vcb.attribute(VcfPeekAf.this.frequencyTag, optFreq.getAsDouble());
				}
			return addFilters3(ctx, vcb, foundLt);
			}

		


		VariantContext addFilters(final VariantContext ctx,final Map<Allele,Double> alt2freq) {			
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			
			final List<Double> af_vector = new ArrayList<>(alt2freq.size());
			boolean foundLt=false;
			for(final Allele alt_allele: ctx.getAlternateAlleles()) {
				final double af;
				if(alt_allele.equals(Allele.SPAN_DEL)) {
					af=0;
					}
				else
					{
					af= alt2freq.getOrDefault(alt_allele, 0.0);
					if(VcfPeekAf.this.af_minimum <= af && af<=VcfPeekAf.this.af_maximum ) {
						foundLt=true;
						}
					}
				af_vector.add(af);
				}
			
			if(!StringUtils.isBlank(VcfPeekAf.this.frequencyTag)) vcb.attribute(VcfPeekAf.this.frequencyTag, af_vector);
			
			return addFilters3(ctx, vcb, foundLt);
			}
		}
	
	/** Use INFO/AC and INFO/AN to get frequency */
	private class InfoAcAnPeeker extends AFPeeker
		{
		@Override
		String getName() { return "ACAN";}
		@Override
		String getDescription() { return "use :  'INFO/AC' divided by 'INFO/AN' to compute the Allele frequency";}
		@Override
		void initialize(final VCFHeader h) {
			if(h.getInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY)==null) {
				throw new IllegalArgumentException("Cannot find INFO="+VCFConstants.ALLELE_COUNT_KEY+" in "+resourceVcfFile);
				}
			if(h.getInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY)==null) {
				throw new IllegalArgumentException("Cannot find INFO="+VCFConstants.ALLELE_NUMBER_KEY+" in "+resourceVcfFile);
				}
			}
		@Override
		final VariantContext sanitize(final VariantContext ctx) {
			final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.noGenotypes();
			if(!peek_variant_id) vcb.noID();
			vcb.unfiltered();
			for(final String key:new HashSet<>(ctx.getAttributes().keySet())) {
				if(key.equals(VCFConstants.ALLELE_COUNT_KEY)) continue;
				if(key.equals(VCFConstants.ALLELE_NUMBER_KEY)) continue;
				vcb.rmAttribute(key);
				}
			return vcb.make();
			}
		
		@Override
		VariantContext applyIgnoringAlt(final VariantContext ctx, final List<VariantContext> overlappers) {
			OptionalDouble optFreq = OptionalDouble.empty();
			for(final VariantContext ctx2:overlappers) {
				if(!ctx2.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
					continue;
					}
				if(!ctx2.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY)) {
					continue;
					}
				final double an = ctx2.getAttributeAsDouble(VCFConstants.ALLELE_NUMBER_KEY,0);
				if(an==0)  {
					continue;
					}
				for(final Double ac:ctx2.getAttributeAsDoubleList(VCFConstants.ALLELE_COUNT_KEY,0.0)) {
					if(ac==null) continue;
					final double af = ac/an;
					if(!optFreq.isPresent()|| optFreq.getAsDouble()>af) {
						optFreq = OptionalDouble.of(af);
						}
					}
				}
			return addFiltersIgnoreAlt(ctx,optFreq);
			}
		
		@Override
		VariantContext applyAlt(final VariantContext ctx, final List<VariantContext> overlappers) {
			final List<Allele> alt_alleles = ctx.getAlternateAlleles();
			final Map<Allele,Double> allele2freq = new HashMap<>(alt_alleles.size());
			
			
			// loop over each ALT from this ctx
			for(int i=0;i< alt_alleles.size();++i) {
				final Allele ctx_alt= alt_alleles.get(i);
				if(ctx_alt.equals(Allele.SPAN_DEL)) {
					continue;
					}
				
				// loop over each variant from database
				for(final VariantContext ctx2:overlappers) {
					if(!ctx2.hasAllele(ctx_alt)) continue;
					
					if(!ctx2.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
						continue;
						}
					if(!ctx2.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY)) {
						continue;
						}
					
					final double an = ctx2.getAttributeAsDouble(VCFConstants.ALLELE_NUMBER_KEY,0);
					if(an==0)  {
						continue;
						}
				
					
					final int ref_idx = ctx2.getAlleleIndex(ctx_alt);
					if(ref_idx<1 /* 0 is ref */)   {
						continue;
						};
					
					final List<Double> ac2 = ctx2.getAttributeAsDoubleList(VCFConstants.ALLELE_COUNT_KEY, 1.0);
					
					if(ref_idx-1>=ac2.size())   {
						continue;
						};
					
					final double alt2freq = ac2.get(ref_idx-1)/an;
					if(allele2freq.getOrDefault(ctx_alt,1.0) > alt2freq) {
						allele2freq.put(ctx_alt, alt2freq);
						}
					}
				}
			return addFilters(ctx,allele2freq);
			}
		}
	
	private abstract class AbstractAFInfoFieldPeeker extends AFPeeker
		{
		protected abstract String getPeekInfoTagName();
		@Override
		void initialize(final VCFHeader h) {
			if(StringUtils.isBlank(getPeekInfoTagName())) {
				throw new IllegalStateException("No INFO tag was defined for extractor "+ this.getName());
				}
			final VCFInfoHeaderLine hl= h.getInfoHeaderLine(getPeekInfoTagName());
			if(hl==null) {
				throw new IllegalArgumentException("Cannot find INFO="+getPeekInfoTagName()+" in "+resourceVcfFile);
				}
			if(!hl.getCountType().equals(VCFHeaderLineCount.A)) {
				LOG.warn("Expected find INFO="+getPeekInfoTagName()+" Count="+VCFHeaderLineCount.A+" but got "+hl.getCountType());
				}
			if(hl.getType().equals(VCFHeaderLineType.Float)) {
				throw new IllegalArgumentException("Expected find INFO="+getPeekInfoTagName()+" Type="+VCFHeaderLineType.Float+" but got "+hl.getType());
				}
			}
		
		@Override
		final VariantContext sanitize(final VariantContext ctx) {
			final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.noGenotypes();
			if(!peek_variant_id) vcb.noID();
			vcb.unfiltered();
			for(final String key:new HashSet<>(ctx.getAttributes().keySet())) {
				if(key.equals(getPeekInfoTagName())) continue;
				vcb.rmAttribute(key);
				}
			return vcb.make();
			}
		
		@Override
		VariantContext applyIgnoringAlt(final VariantContext ctx,final List<VariantContext> overlappers) {
			final OptionalDouble optFreq = overlappers.
					stream().
					filter(ctx2->ctx2.hasAttribute(getPeekInfoTagName())).
					flatMap(ctx2->ctx2.getAttributeAsDoubleList(getPeekInfoTagName(),0.0).stream()).
					mapToDouble(V->V.doubleValue()).
					min();
			
			return addFiltersIgnoreAlt(ctx,optFreq);
			}

		
		@Override
		VariantContext applyAlt(final VariantContext ctx, final List<VariantContext> overlappers) {
			final List<Allele> alt_alleles = ctx.getAlternateAlleles();
			final Map<Allele,Double> allele2freq = new HashMap<>(alt_alleles.size());
			
			
			
			// loop over each ALT from this ctx
			for(int i=0;i< alt_alleles.size();++i) {
				final Allele ctx_alt= alt_alleles.get(i);
				if(ctx_alt.equals(Allele.SPAN_DEL)) {
					continue;
					}
				
				// loop over each variant from database
				for(final VariantContext ctx2:overlappers) {
					if(!ctx2.hasAllele(ctx_alt)) continue;
					
					if(!ctx2.hasAttribute(getPeekInfoTagName())) {
						continue;
						}
					
					
					final int ref_idx = ctx2.getAlleleIndex(ctx_alt);
					if(ref_idx<1 /* 0 is ref */)   {
						continue;
						};
					
					final List<Double> af_array = ctx2.getAttributeAsDoubleList(getPeekInfoTagName(), 1.0);
					
					if(ref_idx-1>=af_array.size())   {
						continue;
						};
					
					final double alt2freq = af_array.get(ref_idx-1);
					if(allele2freq.getOrDefault(ctx_alt,1.0) > alt2freq) {
						allele2freq.put(ctx_alt, alt2freq);
						}
					}
				}
			return addFilters(ctx,allele2freq);
			}
		}
	
	
	private class InfoAfPeeker extends AbstractAFInfoFieldPeeker
		{
		@Override
		String getName() {
			return getPeekInfoTagName();
			}
		@Override
		String getDescription() { return "use  field 'INFO/AF' to extract the Allele frequency";}

		@Override
		protected String getPeekInfoTagName() {
			return VCFConstants.ALLELE_FREQUENCY_KEY;
			}
		}
	
	private class CustomInfoPeeker extends AbstractAFInfoFieldPeeker
		{
		@Override
		String getName() {
			return "Custom";
			}
		@Override
		String getDescription() { return "use ratio INFO field defined with option -P to compute the Allele frequency";}

		@Override
		protected String getPeekInfoTagName() {
			return VcfPeekAf.this.custom_peek_info_name;
			}
		}
	
	private class GtPeeker extends AFPeeker
		{
		@Override
		String getName() { return "GT";}
		@Override
		String getDescription() { return "compute Allele frequency from the Genotypes";}
		
		@Override
		final VariantContext sanitize(final VariantContext ctx) {
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx).unfiltered().attributes(Collections.emptyMap());
			if(!peek_variant_id) vcb.noID();
			return vcb.make();
			}
		
		@Override
		void initialize(final VCFHeader h) {
			if(!h.hasGenotypingData()) {
				throw new IllegalArgumentException("No Genotype in "+resourceVcfFile);
				}
			}
		
		@Override
		VariantContext applyIgnoringAlt(VariantContext ctx, final List<VariantContext> overlappers) {
			OptionalDouble optFreq = OptionalDouble.empty();
			for(final VariantContext ctx2:overlappers) {
				final double an= ctx2.getGenotypes().stream().filter(G->G.isCalled()).mapToInt(G->G.getAlleles().size()).sum();
				if(an==0)  {
					continue;
					}

				for(final Allele alt: ctx2.getAlternateAlleles()) {

					double ac= ctx2.getGenotypes().stream().filter(G->G.isCalled()).flatMap(G->G.getAlleles().stream()).filter(A->A.equals(alt)).count();
					
					
					final double af = ac/an;
					if(!optFreq.isPresent() || optFreq.getAsDouble()>af) optFreq = OptionalDouble.of(af);
					}
				}
			return addFiltersIgnoreAlt(ctx,optFreq);
			}

		
		@Override
		VariantContext applyAlt(final VariantContext ctx, final List<VariantContext> overlappers) {
			final List<Allele> alt_alleles = ctx.getAlternateAlleles();
			final Map<Allele,Double> allele2freq = new HashMap<>(alt_alleles.size());
			// loop over each variant from database
			for(int i=0;i< alt_alleles.size();++i) {
				final Allele ctx_alt= alt_alleles.get(i);
				
				if(ctx_alt.equals(Allele.SPAN_DEL)) {
					continue;
					}
				
				for(final VariantContext ctx2:overlappers) {
				// loop over each ALT from this ctx
					if(!ctx2.hasAllele(ctx_alt)) continue;
					
					
				
					final double an= ctx2.getGenotypes().stream().filter(G->G.isCalled()).mapToInt(G->G.getAlleles().size()).sum();
					if(an==0)  {
						continue;
						}

					final  double ac= ctx2.getGenotypes().stream().filter(G->G.isCalled()).flatMap(G->G.getAlleles().stream()).filter(A->A.equals(ctx_alt)).count();
					
					
					final double alt2freq = ac/an;
					if(allele2freq.getOrDefault(ctx_alt,1.0) > alt2freq) {
						allele2freq.put(ctx_alt, alt2freq);
						}
					}
				}
			return addFilters(ctx,allele2freq);
			}
		}

	private BufferedVCFReader indexedVcfFileReader=null;
	private AFPeeker peeker;
	
	public VcfPeekAf()
		{
		}
	
	
		
	private List<VariantContext> getOverlappingBuffer(
			final String dbContig,
			final VariantContext userCtx
			) {
		try(CloseableIterator<VariantContext> iter = this.indexedVcfFileReader.query(dbContig, userCtx.getStart(), userCtx.getEnd())) {
			return iter.stream().collect(Collectors.toList());
			}
		}
	
	private boolean alleles_match_for_id(final VariantContext userCtx,final VariantContext databaseV) {
		if(userCtx.hasID()) return false;
		if(!databaseV.hasID()) return false;
		if(userCtx.getStart()!=databaseV.getStart()) return false;
		
		final Set<Allele> set1 = new HashSet<>(userCtx.getAlleles());
		set1.remove(Allele.NO_CALL);
		set1.remove(Allele.SPAN_DEL);
		final Set<Allele> set2 = new HashSet<>(databaseV.getAlleles());
		set2.remove(Allele.NO_CALL);
		set2.remove(Allele.SPAN_DEL);
		return set2.containsAll(set1);
		}
	
	@Override
	public int doVcfToVcf(
			final String inputName, 
			final VCFIterator vcfIn,
			final VariantContextWriter out)
		{
		try
			{
			final VCFHeader h = vcfIn.getHeader();
			final SAMSequenceDictionary dictDatabase = this.indexedVcfFileReader.getHeader().getSequenceDictionary();  
			final ContigNameConverter dbCtgConverter = dictDatabase==null || dictDatabase.isEmpty()?
					ContigNameConverter.getIdentity():
					ContigNameConverter.fromOneDictionary(dictDatabase)
					;
			
			final VCFHeader h2 = new VCFHeader(h);
			if(!StringUtils.isBlank(this.frequencyTag)) {
				final String msg="Allele Frequency found in "+resourceVcfFile+" with peeker: "+this.peeker.getName();
				if(this.disable_alt_concordance) {
					h2.addMetaDataLine(new VCFInfoHeaderLine(this.frequencyTag,1,VCFHeaderLineType.Float,"Min "+msg));
					} 
				else {
					h2.addMetaDataLine(new VCFInfoHeaderLine(this.frequencyTag,VCFHeaderLineCount.A,VCFHeaderLineType.Float,msg));
					}
			}
			if(!StringUtils.isBlank(this.filterStr)) {
				h2.addMetaDataLine(new VCFFilterHeaderLine(this.filterStr,"Allele Frequency found in "+resourceVcfFile+" with peeker failing the following state: "+					
						this.af_minimum + " <= "+ this.peeker.getName()+" <= "+this.af_maximum));
			}
			
			
			JVarkitVersion.getInstance().addMetaData(this, h2);
			
			
			
			out.writeHeader(h2);
			final ProgressFactory.Watcher<VariantContext> progress = 
					ProgressFactory.newInstance().
					dictionary(h).
					logger(LOG).
					build()
					;
			while(vcfIn.hasNext())
				{
				final VariantContext ctx=progress.apply(vcfIn.next());
				final String dbContig = dbCtgConverter.apply(ctx.getContig());

				final List<VariantContext> overlappers;
				
				if(StringUtils.isBlank(dbContig))
					{
					overlappers = Collections.emptyList();
					}
				else
					{
					overlappers = this.getOverlappingBuffer(dbContig,ctx);
					}
				
				VariantContext ctx2 = this.peeker.apply(ctx,overlappers);
				if(ctx2==null) continue;
				
				/* peek variant ID */
				if(this.peek_variant_id && !ctx2.hasID() && !overlappers.isEmpty()) {
					final VariantContext ctx2_final = ctx2;
					final String id = overlappers.stream().
						filter(V-> alleles_match_for_id(ctx2_final,V)).
						map(V->V.getID()).		
						findFirst().
						orElse(null);
					if(!StringUtils.isBlank(id)) {
						ctx2 = new VariantContextBuilder(ctx).id(id).make();
						}
					}
				
				out.add(ctx2);
				}
			progress.close();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}

	@Override
	protected int beforeVcf() {
		final List<AFPeeker> all_peekers = new ArrayList<>();
		all_peekers.add(new InfoAcAnPeeker());
		all_peekers.add(new InfoAfPeeker());
		all_peekers.add(new GtPeeker());
		all_peekers.add(new CustomInfoPeeker());
		this.indexedVcfFileReader = null;
		if(this.buffer_size<1) {
			LOG.error("bad buffer-size");
			return -1;
			}
		try
			{
			if(this.list_peekers) {
				try(PrintWriter out=super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				for(final AFPeeker p:all_peekers) out.println(p.getName()+"\n\t"+p.getDescription());
				out.flush();
				}
			System.exit(0);
			}
			
			if(StringUtils.isBlank(this.peekerName)) {
				LOG.error("peeker name is empty");
				return -1;
				}
			
			this.peeker = all_peekers.stream().filter(P->P.getName().equals(this.peekerName)).findFirst().orElse(null);
			if(this.peeker==null) {
				LOG.error("peeker "+this.peekerName+" not found in "+ 
						all_peekers.
							stream().
							map(P->P.getName()).
							collect(Collectors.joining(";")));
				return -1;
				}
						
			final VCFReader reader0 = VCFReaderFactory.makeDefault().open(this.resourceVcfFile,true);
			this.indexedVcfFileReader = new BufferedVCFReader(reader0,this.buffer_size);
			this.peeker.initialize(this.indexedVcfFileReader.getHeader());
			this.indexedVcfFileReader.setSimplifier(peeker::sanitize);
			return 0;
			} 
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	@Override
	protected void afterVcf() {
		CloserUtil.close(this.indexedVcfFileReader);
		this.indexedVcfFileReader=null;
		}
	
	
	public static void main(final String[] args) throws IOException
		{
		new VcfPeekAf().instanceMainWithExit(args);
		}
}
