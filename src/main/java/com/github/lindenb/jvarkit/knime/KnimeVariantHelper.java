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
package com.github.lindenb.jvarkit.knime;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.IndexedVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class KnimeVariantHelper {
	public static final Logger LOG = Logger.build(KnimeVariantHelper.class).make();
	private VepPredictionParser vepPredictionParser=null;
	private final AnnPredictionParser annParser;
	private final Map<String,IndexedBedReader> bedReaders=new HashMap<>();
	private final Map<String,IndexedVcfFileReader> vcfReaders=new HashMap<>();
	
	
	public KnimeVariantHelper() {
		final VCFHeader header=new VCFHeader();
		final VCFInfoHeaderLine info=new VCFInfoHeaderLine(
				AnnPredictionParser.getDefaultTag(),
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				""
				);
		header.addMetaDataLine(info);
		this.annParser=new AnnPredictionParserFactory(header).get();	
		}

	public void dispose()
		{
		for(IndexedBedReader r:this.bedReaders.values()) CloserUtil.close(r);
		for(IndexedVcfFileReader r:this.vcfReaders.values()) CloserUtil.close(r);
		}	
	
	@Override
	protected void finalize() throws Throwable {
		dispose();
		super.finalize();
		}
	
	public IndexedBedReader openBed(final String resourceName,String path) throws IOException {
		failIf(this.bedReaders.containsKey(resourceName), "duplicate resource "+resourceName);
		final IndexedBedReader reader = new IndexedBedReader(path);
		this.bedReaders.put(resourceName, reader);
		return reader;
		}
	public IndexedVcfFileReader openVcf(final String resourceName,String path) throws IOException {
		failIf(this.vcfReaders.containsKey(resourceName), "duplicate resource "+resourceName);
		final IndexedVcfFileReader reader = new IndexedVcfFileReader(path);
		this.vcfReaders.put(resourceName,  reader);
		return reader;
		}
	
	public List<BedLine> getBedLines(final String resourceName,final VariantContext ctx) throws IOException
		{
		IndexedBedReader bedReader = this.bedReaders.get(resourceName);
		if(bedReader==null) {
			LOG.error("No such bedReader "+resourceName);
			return Collections.emptyList();
			}
		return bedReader.getLines(ctx.getContig(),ctx.getStart()-1,ctx.getEnd());
		}
	
	public List<VariantContext> getVcfLines(final String resourceName,final VariantContext ctx) throws IOException
		{
		IndexedVcfFileReader bedReader = this.vcfReaders.get(resourceName);
		if(bedReader==null) {
			LOG.error("No such vcfReader "+resourceName);
			return Collections.emptyList();
			}
		return bedReader.getVariants(ctx.getContig(),ctx.getStart(),ctx.getEnd());
		}

	
	public KnimeVariantHelper initSnpEffParser(final String definition)
		{
		failIf(definition==null || definition.isEmpty(),"SNpeff definition is empty");
		final VCFHeader header=new VCFHeader();
		final VCFInfoHeaderLine info=new VCFInfoHeaderLine(
				VepPredictionParser.getDefaultTag(),
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				definition
				);
		header.addMetaDataLine(info);
		this.vepPredictionParser=new VepPredictionParserFactory(header).get();
		return this;
		}
	
	private static void failIf(boolean testFailed, String msg) {
		if (!testFailed)
			return;
		LOG.error("FAILURE: " + msg);
		throw new RuntimeException("The following test failed: " + msg);
	}

	private static int parseInt(final Object o) {
		if (o == null) {
			LOG.error("null input for integer");
			throw new NumberFormatException("Cannot convert null to integer");
		} else if (o instanceof Integer) {
			return Integer.class.cast(o);
		} else {
			try {
				return Integer.parseInt(String.valueOf(o));
			} catch (NumberFormatException err) {
				LOG.error("Cannot cast \"" + o + "\" to integer", err);
				throw err;
			}
		}
	}

	private static final int[] decodeInts(final String string) {
		List<String> split = ParsingUtils.split(string, ',');
		int[] values = new int[split.size()];
		try {
			for (int i = 0; i < values.length; i++) {
				values[i] = Integer.parseInt(split.get(i));
			}
		} catch (final NumberFormatException e) {
			return null;
		}
		return values;
	}

	public static class VariantBuilder {
		private static class Genot {
			String sample;
			String content;
		}

		private final Pattern commaRegex = Pattern.compile("[,]");
		private final Pattern semiColon = Pattern.compile("[;]");
		private final Pattern colonRegex = Pattern.compile("[:]");
		private String contig = null;
		private String format = null;
		private String filter = null;
		private Integer pos = null;
		private String id = null;
		private Double qual = null;
		private Allele ref = null;
		private final Map<String, Object> attributes = new HashMap<>();
		private final List<Allele> alts = new ArrayList<>();
		private final List<Genot> genotypes = new ArrayList<>();

		public VariantBuilder contig(final String s) {
			failIf(s == null || s.isEmpty(), "contig null or empty");
			this.contig = (s == null ? null : String.valueOf(s));
			return this;
		}
		
		public VariantBuilder contig(final Integer s) {
			return this.contig(s==null?(String)null:String.valueOf(s));
		}

		public VariantBuilder chrom(final String s) {
			return contig(s);
		}

		public VariantBuilder chrom(final Integer s) {
			return contig(s);
		}
		
		public VariantBuilder pos(final Integer s) {
			failIf(s == null || s < 0, "position is null or < 0");
			this.pos = (s == null ? null : parseInt(s));
			return this;
		}

		public VariantBuilder id(final String s) {
			this.id = null;
			if (s == null || s.isEmpty() || s.equals(VCFConstants.EMPTY_ID_FIELD))
				return this;
			this.id = s;
			return this;
		}

		public VariantBuilder filter(final String s) {
			this.filter = null;
			if (s == null || s.isEmpty() || s.equals(VCFConstants.UNFILTERED))
				return this;
			this.filter = s;
			return this;
		}

		public VariantBuilder ref(final String s) {
			failIf(s == null || s.isEmpty(), "Reference null or empty");
			this.ref = Allele.create(s, true);
			return this;
		}

		public VariantBuilder alts(final String s) {
			this.alts.clear();
			if (s == null || s.isEmpty())
				return this;
			if (s.isEmpty())
				return this;
			for (final String as : commaRegex.split(s)) {
				this.alts.add(Allele.create(as, false));
			}
			return this;
		}

		public VariantBuilder qual(final Object o) {
			this.qual = null;
			if (o == null)
				return this;
			final String s = o.toString();
			if (s.isEmpty() || s.equals(".") || s.equals(VCFConstants.MISSING_VALUE_v4))
				return this;
			final Double val = Double.valueOf(s);

			// check to see if they encoded the missing qual score in VCF 3
			// style, with either the -1 or -1.0. check for val < 0 to save some
			// CPU cycles
			if ((val < 0)
					&& (Math.abs(val - VCFConstants.MISSING_QUALITY_v3_DOUBLE) < VCFConstants.VCF_ENCODING_EPSILON)) {
				this.qual = null;
			} else {
				this.qual = val;
			}

			return this;
		}

		public final VariantBuilder info(final String s)
			{	
			return this.attributes(s);
			}
		public VariantBuilder attributes(final String s) {
			this.attributes.clear();
			if (s == null || s.isEmpty())
				return this;
			if (s.isEmpty())
				return this;
			for (final String field : semiColon.split(s)) {
				if (field.isEmpty())
					continue;
				int eq = field.indexOf('=');
				final String key;
				final Object value;
				if (eq == -1) {
					key = field;
					value = Boolean.TRUE;
				} else {
					final List<Object> L = new ArrayList<>();
					key = field.substring(0, eq);
					for (final String v : this.commaRegex.split(field.substring(eq + 1))) {

						Object v2 = null;
						try {
							if (v2 == null)
								v2 = Integer.parseInt(v);
						} catch (NumberFormatException err2) {
							v2 = null;
						}
						try {
							if (v2 == null)
								v2 = Double.parseDouble(v);
						} catch (NumberFormatException err2) {
							v2 = null;
						}

						if (v2 == null)
							v2 = v;
						L.add(v2);
					}
					if (L.size() == 1) {
						value = L.get(0);
					} else if (L.isEmpty()) {
						throw new RuntimeException("Cannot parse INFO: \"" + field + "\"");
					} else {
						value = L;
					}
				}
				this.attributes.put(key, value);
			}
			return this;
		}

		public VariantBuilder format(final String s) {
			failIf(s == null || s.isEmpty(), "FORMAT null or empty");
			this.format = s;
			return this;
		}

		public VariantBuilder genotype(final String sampleName, final String content) {
			failIf(sampleName == null, "Sample Name is null");
			final Genot g = new Genot();
			g.sample = sampleName;
			g.content = content;
			this.genotypes.add(g);
			return this;
		}

		public final VariantContext make() {
			return this.build();
		}
		
		@SuppressWarnings("deprecation")
		public VariantContext build() {
			failIf(this.contig == null, "Contig undefined");
			failIf(this.pos == null, "Position undefined");
			failIf(this.ref == null, "Reference undefined");
			final VariantContextBuilder vcb = new VariantContextBuilder();
			vcb.chr(this.contig);
			vcb.start(this.pos);
			vcb.stop(this.pos + this.ref.length() - 1);
			if (this.id != null)
				vcb.id(this.id);
			final List<Allele> alleles = new ArrayList<>(this.alts.size() + 1);
			alleles.add(this.ref);
			alleles.addAll(this.alts);
			vcb.alleles(alleles);
			vcb.attributes(this.attributes);

			if (this.qual != null) {
				vcb.log10PError(this.qual);
			}

			if (!this.genotypes.isEmpty()) {
				final List<Genotype> gtList = new ArrayList<>();
				failIf(this.format == null || this.format.isEmpty(), "Genotypes defined but not FORMAT");
				final String formatTokens[] = colonRegex.split(this.format);
				for (final Genot gt : this.genotypes) {
					if (gt.content == null || gt.content.equals(".") || gt.content.isEmpty()) {
						gtList.add(GenotypeBuilder.createMissing(gt.sample, 2));
					} else {
						final String gtTokens[] = colonRegex.split(gt.content);
						final GenotypeBuilder gb = new GenotypeBuilder(gt.sample);

						for (int i = 0; i < gtTokens.length && i < formatTokens.length; ++i) {
							final String gtKey = formatTokens[i];
							if (gtTokens[i].equals(".") || gtTokens[i].isEmpty()) {
								continue;
							} else if (gtKey.equals(VCFConstants.GENOTYPE_KEY)) {
								final StringTokenizer st = new StringTokenizer(gtTokens[i],
										VCFConstants.PHASING_TOKENS);
								List<Allele> GTAlleles = new ArrayList<Allele>(st.countTokens());
								while (st.hasMoreTokens()) {
									final String indexStr = st.nextToken();

									if (indexStr.equals(VCFConstants.EMPTY_ALLELE)) {
										GTAlleles.add(Allele.NO_CALL);
									} else {
										int j = parseInt(indexStr);

										if (j >= 0 && j < alleles.size()) {
											GTAlleles.add(alleles.get(j));
										} else {
											throw new RuntimeException(
													"Illegal Allele index " + indexStr + " for " + alleles);
										}
									}
								}
								gb.alleles(GTAlleles);
							} else if (gtKey.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS)) {
								gb.AD(decodeInts(gtTokens[i]));
							} else if (gtKey.equals(VCFConstants.GENOTYPE_PL_KEY)) {
								gb.PL(decodeInts(gtTokens[i]));
							} else if (gtKey.equals(VCFConstants.GENOTYPE_LIKELIHOODS_KEY)) {
								gb.PL(GenotypeLikelihoods.fromGLField(gtTokens[i]).getAsPLs());
							} else if (gtKey.equals(VCFConstants.DEPTH_KEY)) {
								gb.DP(Integer.valueOf(gtTokens[i]));
							} else {
								gb.attribute(gtKey, gtTokens[i]);
							}
						}
						gtList.add(gb.make());
					}
				}
				vcb.genotypes(gtList);
			}

			if (this.filter != null) {
				if (this.filter.equals(VCFConstants.PASSES_FILTERS_v4)
						|| this.filter.equals(VCFConstants.PASSES_FILTERS_v3)) {
					vcb.passFilters();
				} else {
					vcb.filters(Arrays.stream(this.commaRegex.split(this.filter)).collect(Collectors.toSet()));
				}
			}

			try {
				return vcb.make();
			} catch (final Throwable err) {
				LOG.error("Cannot convert to variant", err);
				throw err;
			}
		}
	}

	public VariantBuilder build() {
		return new VariantBuilder();
	}
	
	public SequenceOntologyTree getSequenceOntologyTree() {
		return  SequenceOntologyTree.getInstance();
	}

	public List<AnnPredictionParser.AnnPrediction> getAnnPredictions(final VariantContext ctx) {
		return this.annParser.getPredictions(ctx);
		}
	
	public List<VepPredictionParser.VepPrediction> getVepPredictions(final VariantContext ctx) {
		if(this.vepPredictionParser==null) return Collections.emptyList();
		return this.vepPredictionParser.getPredictions(ctx);
		}

	/** return true if variant has any prediction with a SO term (or its children) with this label */
	public boolean hasSequenceOntologyLabel(final VariantContext ctx,final String lbl)
		{
		if(lbl==null) return false;
		final SequenceOntologyTree.Term t= this.getSequenceOntologyTree().getTermByLabel(lbl);
		if(t==null) LOG.warning("don't know SO.label "+lbl);
		return hasSequenceOntologyTerm(ctx,t);
		}
	/** return true if variant has any prediction with a SO term (or its children) with this accession */
	public boolean hasSequenceOntologyAccession(final VariantContext ctx,final String acn)
		{
		if(acn==null) return false;
		final SequenceOntologyTree.Term t= this.getSequenceOntologyTree().getTermByAcn(acn);
		if(t==null) LOG.warning("don't know SO.acn "+acn);
		return hasSequenceOntologyTerm(ctx,t);
		}

	/** return true if variant has any prediction with a SO term (or its children) */
	public boolean hasSequenceOntologyTerm(final VariantContext ctx,final SequenceOntologyTree.Term t)
		{
		if(t==null) return false;
		final Set<SequenceOntologyTree.Term> children=t.getAllDescendants();
		for(AnnPredictionParser.AnnPrediction a: getAnnPredictions(ctx)) {
			if(!Collections.disjoint(a.getSOTerms(),children)) return true;
			}
		for(VepPredictionParser.VepPrediction a: getVepPredictions(ctx)) {
			if(!Collections.disjoint(a.getSOTerms(),children)) return true;
			}
		return false;
		}
	
	public boolean isMendelianIncompatibility(final Genotype child,final Genotype parent)
		{
		if(child==null || parent==null) return false;
		if(child.isNoCall() || parent.isNoCall()) return false;
		if(child.getPloidy()!=2 || parent.getPloidy()!=2) return false;
		for(final Allele childAllele:child.getAlleles())
			{
			if(parent.getAlleles().contains(childAllele)) return false;
			}
		
		return true;
		}
	public boolean isMendelianIncompatibility(final Genotype child,final Genotype father,final Genotype mother)
		{
		if(child==null || child.isNoCall()) return false;
		if(father==null || father.isNoCall()) {
			return this.isMendelianIncompatibility(child,mother);
			}
		if(mother==null || mother.isNoCall()) {
			return this.isMendelianIncompatibility(child,father);
			}
		final Allele alleles[]=new Allele[2];
		for(final Allele af:father.getAlleles())
		{
			alleles[0]=af;
			for(final Allele am:mother.getAlleles())
			{
				alleles[1]=am;
				final Genotype sim = new GenotypeBuilder(child.getSampleName()).alleles(Arrays.asList(alleles)).make();
				if(child.sameGenotype(sim, true)) return false;
			}	
		}
		
		return true;
		}

	
}
