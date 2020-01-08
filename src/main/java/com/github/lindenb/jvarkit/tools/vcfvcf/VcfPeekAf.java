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
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.variant.vcf.VCFIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
/**
BEGIN_DOC

END_DOC

 */
@Program(name="vcfpeekaf",
		description="Peek the AF from another VCF",
		keywords={"vcf","annotation","af"},
		modificationDate="20191202",
		creationDate="20191202"
		)
public class VcfPeekAf extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfPeekAf.class).make();

	
	@Parameter(names={"-F","--database","--tabix","--resource"},description="An indexed VCF file. Source of the annotations",required=true)
	private Path resourceVcfFile = null;	
	@Parameter(names={"-T","--tag"},description="INFO tag to put found frequency. empty: no extra tag.")
	private String frequencyTag = "";
	@Parameter(names={"-f","--filter"},description="soft FILTER the variant of this data if AF is not found or it greater > threshold. If empty, just DISCARD the variant")
	private String filterStr = "";
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-b","--buffer-size"},converter=DistanceParser.StringConverter.class, description="buffer size (in bp). We don't do a random access for each variant. Instead of this, load all the variants in a defined window. "+DistanceParser.OPT_DESCRIPTION,splitter=NoSplitter.class)
	private int buffer_size = 100_000;
	@Parameter(names={"-l","--list"},description="List available AF peekers and exit.",help=true)
	private boolean list_peekers = false;
	@Parameter(names={"-p","--peeker"},description="Peeker name",required=true)
	private String peekerName = null;
	@Parameter(names={"-t","--treshold"},description="AF treshold. Variant is accepted is computed AF <= treshold.",required=true)
	private double af_treshold = 1.0;
	@Parameter(names={"--params"},description="Extra parameters for the AF peeker",hidden=true)
	private List<String> peeker_parameters = new ArrayList<>();
	@Parameter(names={"--no-alt"},description="Do not look at the alternate alleles concordance")
	private boolean disable_alt_concordance = false;


	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	
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
		abstract VariantContext applyAlt(final VariantContext ctx,final List<VariantContext> overlappers);
		abstract VariantContext applyIgnoringAlt(final VariantContext ctx,final List<VariantContext> overlappers);
		
		final VariantContext apply(final VariantContext ctx,final List<VariantContext> overlappers) {
			return disable_alt_concordance ? applyIgnoringAlt(ctx,overlappers):applyAlt(ctx, overlappers);
			}

		
		abstract VariantContext sanitize(final VariantContext ctx);
		
		VariantContext addFilters3(final VariantContext ctx,final VariantContextBuilder vcb,boolean accept) {
			final Set<String> old_filters = ctx.getFilters();
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

		
		VariantContext addFiltersIgnoreAlt(final VariantContext ctx,final Double minFreq) {
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			final boolean foundLt=minFreq==null || minFreq.doubleValue()< VcfPeekAf.this.af_treshold;
			if(minFreq!=null && !StringUtils.isBlank(VcfPeekAf.this.frequencyTag)) vcb.attribute(VcfPeekAf.this.frequencyTag, minFreq);
			return addFilters3(ctx, vcb, foundLt);
			}

		


		VariantContext addFilters(final VariantContext ctx,final Map<Allele,Double> alt2freq) {			
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			
			final List<Double> af_vector = new ArrayList<>();
			boolean foundLt=false;
			for(final Allele alt_allele: ctx.getAlternateAlleles()) {
				final double af;
				if(alt_allele.equals(Allele.SPAN_DEL)) {
					af=0;
					}
				else
					{
					af= alt2freq.getOrDefault(alt_allele, 0.0);
					if(af<=VcfPeekAf.this.af_treshold ) {
						foundLt=true;
						}
					}
				af_vector.add(af);
				}
			
			if(!StringUtils.isBlank(VcfPeekAf.this.frequencyTag)) vcb.attribute(VcfPeekAf.this.frequencyTag, af_vector);
			
			return addFilters3(ctx, vcb, foundLt);
			}
		}
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
		VariantContext sanitize(final VariantContext ctx) {
			final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.noGenotypes();
			vcb.noID();
			for(final String key:new HashSet<>(ctx.getAttributes().keySet())) {
				if(key.equals(VCFConstants.ALLELE_COUNT_KEY)) continue;
				if(key.equals(VCFConstants.ALLELE_NUMBER_KEY)) continue;
				vcb.rmAttribute(key);
				}
			return vcb.make();
			}
		
		@Override
		VariantContext applyIgnoringAlt(VariantContext ctx, List<VariantContext> overlappers) {
			Double minFreq=null;
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
				for(final Double af:ctx2.getAttributeAsDoubleList(VCFConstants.ALLELE_COUNT_KEY,0.0)) {
					if(minFreq==null || minFreq.doubleValue()>af) minFreq = af;
					}
				}
			return addFiltersIgnoreAlt(ctx,minFreq);
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
	private class InfoAfPeeker extends AFPeeker
		{
		@Override
		String getName() { return "AF";}
		@Override
		String getDescription() { return "use ratio INFO/AF to compute the Allele frequency";}
		@Override
		void initialize(final VCFHeader h) {
			if(h.getInfoHeaderLine(VCFConstants.ALLELE_FREQUENCY_KEY)==null) {
				throw new IllegalArgumentException("Cannot find INFO="+VCFConstants.ALLELE_FREQUENCY_KEY+" in "+resourceVcfFile);
				}
			}
		@Override
		VariantContext sanitize(final VariantContext ctx) {
			final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			vcb.noGenotypes();
			vcb.noID();
			for(final String key:new HashSet<>(ctx.getAttributes().keySet())) {
				if(key.equals(VCFConstants.ALLELE_FREQUENCY_KEY)) continue;
				vcb.rmAttribute(key);
				}
			return vcb.make();
			}
		
		@Override
		VariantContext applyIgnoringAlt(VariantContext ctx, List<VariantContext> overlappers) {
			Double minFreq=null;
			for(final VariantContext ctx2:overlappers) {
				if(!ctx2.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
					continue;
					}
				
				for(final Double af:ctx2.getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY,0.0)) {
					if(minFreq==null || minFreq.doubleValue()>af) minFreq = af;
					}
				}
			return addFiltersIgnoreAlt(ctx,minFreq);
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
					
					if(!ctx2.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
						continue;
						}
					
				
					
					final int ref_idx = ctx2.getAlleleIndex(ctx_alt);
					if(ref_idx<1 /* 0 is ref */)   {
						continue;
						};
					
					final List<Double> af_array = ctx2.getAttributeAsDoubleList(VCFConstants.ALLELE_FREQUENCY_KEY, 1.0);
					
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
	private class GtPeeker extends AFPeeker
		{
		@Override
		String getName() { return "GT";}
		@Override
		String getDescription() { return "compute Allele frequency from the Genotypes";}
		
		@Override
		VariantContext sanitize(VariantContext ctx) {
			return new VariantContextBuilder(ctx).noID().attributes(Collections.emptyMap()).make();
			}
		
		@Override
		void initialize(final VCFHeader h) {
			if(!h.hasGenotypingData()) {
				throw new IllegalArgumentException("No Genotype in "+resourceVcfFile);
				}
			}
		
		@Override
		VariantContext applyIgnoringAlt(VariantContext ctx, List<VariantContext> overlappers) {
			Double minFreq=null;
			for(final VariantContext ctx2:overlappers) {
				double an= ctx2.getGenotypes().stream().filter(G->G.isCalled()).mapToInt(G->G.getAlleles().size()).sum();
				if(an==0)  {
					continue;
					}

				for(final Allele alt: ctx2.getAlternateAlleles()) {

					double ac= ctx2.getGenotypes().stream().filter(G->G.isCalled()).flatMap(G->G.getAlleles().stream()).filter(A->A.equals(alt)).count();
					
					
					final double af = ac/an;
					if(minFreq==null || minFreq.doubleValue()>af) minFreq = af;
					}
				}
			return addFiltersIgnoreAlt(ctx,minFreq);
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

	private VCFFileReader indexedVcfFileReader=null;
	private final List<VariantContext> buffer = new ArrayList<>();
	private Locatable last_buffer_interval = null;
	private AFPeeker peeker;
	
	public VcfPeekAf()
		{
		}
	
	
		
	private List<VariantContext> getOverlappingBuffer(
			final String dbContig,
			final VariantContext userCtx
			) {
		final int start = userCtx.getStart();
		final int end = userCtx.getEnd();
		if(	!(
			this.last_buffer_interval!=null &&
			this.last_buffer_interval.getContig().equals(dbContig) &&
			this.last_buffer_interval.getStart() < start && 
			end < this.last_buffer_interval.getEnd()
			))
			{
			this.buffer.clear();
			
			this.last_buffer_interval = new SimpleInterval(
					dbContig,
					Math.max(0,start-1),
					(end+1+this.buffer_size)
					);
			
			try( CloseableIterator<VariantContext> t = this.indexedVcfFileReader.query(
					dbContig,
					Math.max(0,start-1),
					(end+1+this.buffer_size)
					)) {
				while(t.hasNext())
					{
					this.buffer.add(this.peeker.sanitize(t.next()));
					}
				}
			}
		return this.buffer.stream().
				filter(V->V.getContig().equals(dbContig) && V.getStart()==start && V.getReference().equals(userCtx.getReference())).
				collect(Collectors.toList());
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
			final SAMSequenceDictionary dictDatabase = this.indexedVcfFileReader.getFileHeader().getSequenceDictionary();  
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
				h2.addMetaDataLine(new VCFFilterHeaderLine(this.filterStr,"Allele Frequency found in "+resourceVcfFile+" with peeker: "+this.peeker.getName()+" is > "+this.af_treshold));
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
				
				final VariantContext ctx2 = this.peeker.apply(ctx,overlappers);
				if(ctx2==null) continue;
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
	public int doWork(final List<String> args) {
		final List<AFPeeker> all_peekers = new ArrayList<>();
		all_peekers.add(new InfoAcAnPeeker());
		all_peekers.add(new InfoAfPeeker());
		all_peekers.add(new GtPeeker());
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
						
			this.indexedVcfFileReader = new VCFFileReader(this.resourceVcfFile,true);
			this.peeker.initialize(this.indexedVcfFileReader.getFileHeader());			
			return doVcfToVcfPath(args,this.writingVariantsDelegate, this.outputFile);
			} 
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedVcfFileReader);
			this.indexedVcfFileReader=null;
			}
		}
	
	
	public static void main(final String[] args) throws IOException
		{
		new VcfPeekAf().instanceMainWithExit(args);
		}
}
