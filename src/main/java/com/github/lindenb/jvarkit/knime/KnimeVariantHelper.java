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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.knime;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.misc.VcfToTable;
import com.github.lindenb.jvarkit.tools.vcfvcf.VcfPeekVcf;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.IndexedVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.IndexFactory.IndexType;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

```
 try {
	final com.github.lindenb.jvarkit.knime.KnimeVariantHelper helper = new com.github.lindenb.jvarkit.knime.KnimeVariantHelper();
	
	htsjdk.variant.variantcontext.VariantContext ctx = helper.build().
	   contig($#CHROM$).
	   pos($POS$).
	   id($ID$).
	   ref($REF$).
	   alts($ALT$).
	   filter($FILTER$).
	   info($INFO$).
	   format($FORMAT$).
	   genotype("CD12100",$CD12100$).
	   genotype("CD12218",$CD12218$).
	   build()
	   ;
	
	helper.initSnpEffParser("Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SYMBOL_SOURCE|HGNC_ID|RefSeq|SIFT|PolyPhen");
	
	
	
	return  ctx.hasID()==false  && ctx.isIndel() &&
	        helper.hasSequenceOntologyLabel(ctx,"protein_altering_variant")  && 
	        ctx.getAlternateAlleles().size()==1 && 
			ctx.getAttributeAsDouble("AF",1.0) > 0.0001  &&
	        !ctx.getGenotype("CD12100").sameGenotype( ctx.getGenotype("CD12218") )
	    ;
	}
catch(Throwable err)
	{
	System.err.println("################### ERROR with "+ $#CHROM$ +" "+$POS$);
	err.printStackTrace();
	return false;
	}
	 
```

END_DOC

 */
import htsjdk.variant.vcf.VCFHeaderLine;
@Program(name="knimehelper",description="library for using htsjdk into knime",
		keywords={"knime","vcf"}
		)	
public class KnimeVariantHelper extends VcfTools {
	public static final Logger LOG = Logger.build(KnimeVariantHelper.class).make();
	private final Map<String,IndexedBedReader> bedReaders=new HashMap<>();
	private final Map<String,IndexedVcfFileReader> vcfReaders=new HashMap<>();
	private enum ForceSuffix { No, ForceTabix,ForceTribble};
	private ForceSuffix forceSuffix=ForceSuffix.No;
	private final int SUFFIX_LENGTH=20;
	private File workfingDirectory=null;
	private String fileSuffix=null;
	private String filePrefix=null;
	private final Set<VCFHeaderLine> extraVcfHeaderLines = new HashSet<>();
	/** number of variant printed the last time we called  processVcfMulti */
	private int lastVariantCount=0;
	
	public final VariantContextWriterBuilder variantContextWriterBuilder = 
			new VariantContextWriterBuilder().
				setCreateMD5(false);

	@SuppressWarnings("unused")
	private static VcfPeekVcf __fooljavac1 = null;
	@SuppressWarnings("unused")
	private static VcfToTable __fooljavac2 = null;
	
	public KnimeVariantHelper() {
		
		}
	public KnimeVariantHelper(final VCFHeader header) {
		super(header);
		}

	
	/** called at end to dispose things */
	public void dispose()
		{
		for(final IndexedBedReader r:this.bedReaders.values()) CloserUtil.close(r);
		for(final IndexedVcfFileReader r:this.vcfReaders.values()) CloserUtil.close(r);
		}	
	
	/** get current Logger */
	public Logger getLogger() {
		return LOG;
	}
	
	@Override
	protected void finalize() throws Throwable {
		dispose();
		super.finalize();
		}
	
	
	/** expose the extra vcf header lines, give the user a chance to add some new header line before the next VariantContextWriter.writeHeader */
	public Set<VCFHeaderLine> getExtraVcfHeaderLines()
		{
		return this.extraVcfHeaderLines;
		}
	
	
	public KnimeVariantHelper addVcfHeaderLines(final VCFHeaderLine newheader)
		{
		if(newheader==null) 
			{
			final String msg = "Cannot add newheader==null";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		getExtraVcfHeaderLines().add(newheader);
		return this;
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
				this.qual = val/-10.0;
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
	
	public static class Encoder
		{
		private final VariantContext context;
		private final List<String> tokens;
		Encoder(final VariantContext ctx) {
			this.context = ctx;
			final VCFHeader dummy=new VCFHeader();
			final VCFEncoder encoder = new VCFEncoder(dummy, true,false);
			final String line=encoder.encode(context);
			this.tokens = Arrays.asList(line.split("[\t]"));
			}
		public String get(int index) { return this.tokens.get(index);}
		public String getContig() { return context.getContig();}
		public int getPos() { return context.getStart();}
		public String getId() { return get(2);}
		public String getRef() { return get(3);}
		public String getAlt() { return get(4);}
		public String getQual() { return get(5);}
		public String getFilter() { return get(6);}
		public String getInfo() { return get(7);}
		public String getFormat() { return get(8);}
		public String getGenotype(int i) { return get(9+i);}
		public String getGenotype(final String sampleName) { 
			for(int i=0;i< this.context.getNSamples();++i)
				{
				if(this.context.getGenotype(i).getSampleName().equals(sampleName))
					{
					return getGenotype(i);
					}
				}
			throw new IllegalArgumentException("no such sample :"+sampleName);
			}
		@Override
		public String toString() {
			return String.join("\t", this.tokens);
			}
		}
	
	public Encoder encoder(final VariantContext ctx) {
		return new Encoder(ctx);
	}
	
	/** build a dummy VCF header from Variant Context 
	private VCFHeader makeDummyHeader(final VariantContext ctx) {
		final List<String> sampleNames = new ArrayList<>(ctx.getNSamples());
		final Set<VCFHeaderLine> metaData = new HashSet<>();
		for(int i=0;i< ctx.getNSamples();++i) sampleNames.add(ctx.getGenotype(i).getSampleName());
		for(final String filter:ctx.getFilters())
			{
			metaData.add(new VCFFilterHeaderLine(filter));
			}
		final Set<String> formatKey = ctx.getGenotypes().stream().flatMap(G->G.getExtendedAttributes().keySet().stream()).collect(Collectors.toSet());
		
		return new VCFHeader(metaData,sampleNames);
	}*/

/** index a vcf file is needed */
public void indexVcfFile(final String file) throws IOException{
	indexVcfFile(new File(file));
	}	

/** index a vcf file is needed */
public void indexVcfFile(final File file) throws IOException{
	IOUtil.assertFileIsReadable(file);
	if(file.getName().endsWith(".vcf.gz"))
		{
		LOG.info("writing tabix index for "+file);
		final File output=new File(file.getAbsolutePath()+TabixUtils.STANDARD_INDEX_EXTENSION);
		try
			{
			if(output.exists())
    			{
    			getLogger().info("Tabix index "+output+" already exists.");
    			return;
    			}
			final TabixIndex index=IndexFactory.createTabixIndex(file,new VCFCodec(),(SAMSequenceDictionary)null);
			index.write(output);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			throw new IOException(err);
			}
		}
	else if(file.getName().endsWith(".vcf"))
		{
		LOG.info("writing tribble index for "+file);
		final File output=new File(file.getAbsolutePath()+Tribble.STANDARD_INDEX_EXTENSION);
		try
			{
			if(output.exists())
    			{
				getLogger().info("Tribble index "+output+" already exists.");
    			}
			final Index index=IndexFactory.createIndex(file,new VCFCodec(),IndexType.LINEAR);
			index.writeBasedOnFeatureFile(file);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			throw new IOException(err);
			}
		}
	else
		{
		throw new IOException("Cannot index VCF file "+file);
		}
	}


/** wrapper around Variant , contains header, tools and variant*/
public static abstract class VariantData
	{
	private final VariantContext ctx;
	VariantData(final VariantContext ctx){
		this.ctx = ctx;
		}
	public abstract VCFHeader getVcfHeader();
	public abstract VcfTools getVcfTools();
	public VariantContext getVariant() { return this.ctx;}
	}

public Stream<VariantContext> forEachVariants(final String vcfFile) throws IOException
	{
	final File file = new File(vcfFile);
	IOUtil.assertFileIsReadable(file);
	final VCFFileReader r= new VCFFileReader(file, false);
	final VCFHeader header = r.getFileHeader();
	this.init(header);
	final CloseableIterator<VariantContext> iter = r.iterator();
	final Iterable<VariantContext> iterable = () -> iter;
	return StreamSupport.stream(iterable.spliterator(),false).onClose(()->{
		CloserUtil.close(iter);
		CloserUtil.close(r);
		});
	}

public Stream<VariantData> forEachVariantData(final String vcfFile) throws IOException
	{
	final File file = new File(vcfFile);
	IOUtil.assertFileIsReadable(file);
	final VCFFileReader r= new VCFFileReader(file, false);
	final VCFHeader header = r.getFileHeader();
	this.init(header);
	final VcfTools vcfTools = new VcfTools(header);
	final CloseableIterator<VariantContext> iter = r.iterator();
	final Iterable<VariantContext> iterable = () -> iter;
	return StreamSupport.stream(iterable.spliterator(),false).onClose(()->{
		CloserUtil.close(iter);
		CloserUtil.close(r);
		}).map(V->new VariantData(V){
			@Override public  VCFHeader getVcfHeader() {
				return header;
				}
			@Override public  VcfTools getVcfTools(){
				return vcfTools;
				}
		});
	}


public Predicate<VariantContext> parseVariantIntervalFilters(final String...array)
	{
	final IntervalParser parser=new IntervalParser();
	parser.setRaiseExceptionOnError(true);
	Predicate<VariantContext> filter = V -> false;
	for(final String str:array)	
		{
		final Interval interval = parser.parse(str);
		filter = filter.or(V -> (V.getContig().equals(interval.getContig()) && !(interval.getEnd()<V.getStart() || V.getEnd()<interval.getStart())));
		}
	return filter;
	}

public IntervalTreeMap<Boolean> parseBedAsBooleanIntervalTreeMap(final String bedUri) throws IOException
	{
	if(bedUri==null || bedUri.isEmpty()) throw new IllegalArgumentException("bad bed uri");
	BufferedReader r=null;
	final IntervalTreeMap<Boolean> treeMap = new IntervalTreeMap<>();
	try
		{
		r = IOUtils.openURIForBufferedReading(bedUri);
		final BedLineCodec codec=new BedLineCodec();
		r.lines().forEach(L->{
			final BedLine bedline = codec.decode(L);
			if(bedline==null) {
				LOG.warn("Ignoring line in BED (doesn't look like a BED line):"+L);
				return;
			}
			treeMap.put(bedline.toInterval(),Boolean.TRUE);
			});
		
		return treeMap;
		}
	finally
		{
		CloserUtil.close(r);
		}
	
}






	public KnimeVariantHelper prefix(final String prefix)
		{
		if(prefix!=null && !prefix.isEmpty())
			{
			this.filePrefix = prefix;
			if(!this.filePrefix.endsWith("."))
				{
				this.filePrefix+=".";
				}
			}
		return this;
		}
	
	public KnimeVariantHelper compressed(boolean bgzip)
		{
		this.forceSuffix = (bgzip?ForceSuffix.ForceTabix:ForceSuffix.ForceTribble);
		return this;
		}
	
	/** alias of tmpDir */
	public final KnimeVariantHelper workDir(final String tmpDir)
		{
		return tmpDir(tmpDir);
		}
	/** set the directory where file will be written */
	public KnimeVariantHelper tmpDir(final String tmpDir)
		{
		if(tmpDir!=null && !tmpDir.isEmpty())
			{
			this.workfingDirectory = new File(tmpDir);
			}
		return this;
		}
	
	public KnimeVariantHelper fileId(final String fileId)
		{
		if(fileId!=null)
			{
			this.fileSuffix = fileId;
			while(this.fileSuffix.length()< SUFFIX_LENGTH)
				{
				this.fileSuffix = "0"+this.fileSuffix;
				}
			}
		return this;
		}

	/** creates a generic output file from this input */
	public String createOutputFile(final String pathname,String extension)
		{
		if(pathname==null || pathname.isEmpty())
			{
			final String msg="User Error: pathname was not specified";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		

		if(this.fileSuffix==null || this.fileSuffix.isEmpty())
			{
			final String msg="User Error: file suffix was not specified";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		if(this.fileSuffix.length()!=SUFFIX_LENGTH)
			{
			final String msg="User Error: file suffix (\""+this.fileSuffix+"\") must have length="+this.SUFFIX_LENGTH+". But got "+this.fileSuffix.length()+".";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		if(!this.fileSuffix.matches("[0-9A-Za-z_]{"+SUFFIX_LENGTH+"}"))
			{
			final String msg="User Error: file suffix (\""+this.fileSuffix+"\") must match the regular expression \"[0-9A-Za-z_]{"+SUFFIX_LENGTH+"}\"  .";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		if(this.workfingDirectory==null)
			{
			final String msg="Working directory was not specified.";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		IOUtil.assertDirectoryIsReadable(this.workfingDirectory);
		IOUtil.assertDirectoryIsWritable(this.workfingDirectory);

		final File inVcfFile=new File(pathname);
		String filename= inVcfFile.getName();
		for(;;)
			{
			if(filename.endsWith(".gz"))
				{
				filename = filename.substring(0,filename.length()-3);
				}
			else if(filename.endsWith(".vcf") || 
					filename.endsWith(".txt") || 
					filename.endsWith(".tsv") || 
					filename.endsWith(".csv") || 
					filename.endsWith(".xsl"))
				{
				filename = filename.substring(0,filename.length()-4);
				}
			else
				{
				break;
				}
			}
		if(filename.matches(".*\\.[A-Za-z0-9_]{"+SUFFIX_LENGTH+"}$"))
			{
			filename = filename.substring(0,filename.length()-(SUFFIX_LENGTH+1));
			}
		if((this.filePrefix!=null && !this.filePrefix.isEmpty()) && !filename.startsWith(this.filePrefix))
			{
			filename= this.filePrefix+filename;
			}
		
		if(!extension.startsWith(".")) extension="."+extension;
		filename += "."+this.fileSuffix+ extension;
		
		filename = new File(this.workfingDirectory,filename).getPath();
		
		return filename;
		}
	
	/** 
	 * filter the VCF file,
	 * 
	 * @param vcfIn input file name
	 * @param fun functional
	 * @return the output file name
	 * @throws IOException
	 */
	public String filterVcf(final String vcfIn,final Predicate<VariantContext> filter)
		throws IOException
		{
		if(filter==null )
			{
			final String msg="User Error: predicate was not specified";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		return processVcf(vcfIn,V->filter.test(V)?V:null);
		}
	/** 
	 * process the VCF file,
	 * 
	 * @param vcfIn input file name
	 * @param fun functional
	 * @return the output file name
	 * @throws IOException
	 */
	public String processVcf(final String vcfIn,final Function<VariantContext,VariantContext> fun)
			throws IOException
		{
		if(fun==null )
			{
			final String msg="User Error: function was not specified";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		return processVcfMulti(vcfIn,V->{
			final VariantContext v2 = fun.apply(V);
			return v2==null?
					Collections.emptyList():
					Collections.singletonList(v2);}
				);
		}
	
	/** 
	 * process the VCF file,
	 * 
	 * @param vcfIn input file name
	 * @param fun functional
	 * @return the output file name
	 * @throws IOException
	 */
	public String processVcfMulti(final String vcfIn,
			final Function<VariantContext,List<VariantContext>> fun)
		throws IOException
		{
		this.lastVariantCount=0;
		if(vcfIn==null)
			{
			final String msg="Vcf Input URI/FIle is null.";
			LOG.error(msg);
			throw new IllegalArgumentException(msg);
			}
		
		
		File outVcfFile = null;
		File outVcfIndexFile = null;
		final File STOP_FILE= new File(this.workfingDirectory,"STOP");
		if(STOP_FILE.exists())
			{
			final String msg="There is a stop file in "+STOP_FILE;
			LOG.error(msg);
			throw new IOException(msg);
			}
		boolean fail_flag=false;
		VCFIterator iter=null;
		VariantContextWriter variantContextWriter = null;
		try
			{
			IOUtil.assertDirectoryIsReadable(this.workfingDirectory);
			IOUtil.assertDirectoryIsWritable(this.workfingDirectory);
						
			if(!IOUtil.isUrl(vcfIn)) {
				IOUtil.assertFileIsReadable(new File(vcfIn));
				}
			

			
			final String extension;
			if( this.forceSuffix.equals(ForceSuffix.ForceTabix))
				{
				extension=".vcf.gz";
				}
			else if( this.forceSuffix.equals(ForceSuffix.ForceTribble))
				{
				extension=".vcf";
				}
			else if(vcfIn.endsWith(".gz"))
				{
				extension=".vcf.gz";
				}
			else
				{
				extension=".vcf";
				}
			final String filename  = this.createOutputFile(vcfIn, extension);
			
			final String indexFilename;
			if(extension.endsWith(".gz"))
				{
				indexFilename = filename + Tribble.STANDARD_INDEX_EXTENSION;
				}
			else
				{
				indexFilename = filename + TabixUtils.STANDARD_INDEX_EXTENSION;
				}
			
			outVcfFile = new File(filename);
			outVcfIndexFile = new File(indexFilename);
			LOG.info("opening "+vcfIn);
			iter = VCFUtils.createVCFIterator(vcfIn);
			
			
			
			
			super.init( iter.getHeader());
			
			final VCFHeader vcfHeader2 ;
			if(this.getExtraVcfHeaderLines().isEmpty())
				{
				vcfHeader2 =  iter.getHeader();
				}
			else
				{
				vcfHeader2  =new VCFHeader(iter.getHeader());
				for(final VCFHeaderLine extra: this.getExtraVcfHeaderLines())
					{
					vcfHeader2.addMetaDataLine(extra);
					}
				// clear vcf header line now they 've been added to the header.
				this.getExtraVcfHeaderLines().clear();
				}
				
			
			final SAMSequenceDictionary dict = this.getHeader().getSequenceDictionary();
			if(dict==null)
				{
				final String msg="There is no dictionary (##contig lines) in "+vcfIn+" but they are required.";
				LOG.error(msg);
				throw new IllegalArgumentException(msg);
				}
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(dict);
			progress.setLogPrefix(this.filePrefix);
			
			LOG.info("writing "+outVcfFile+". Emergency stop file is "+STOP_FILE);
			
			
			variantContextWriter = this.variantContextWriterBuilder.
				setOutputFile(outVcfFile).
				setReferenceDictionary(dict).
				build()
				;
			long lastTick = System.currentTimeMillis();
			variantContextWriter.writeHeader(vcfHeader2);
			while(iter.hasNext())
				{
				final VariantContext ctx = progress.watch(iter.next());
				final List<VariantContext> array = fun.apply(ctx);
				if(array!=null)
					{
					for(final VariantContext ctx2: array)
						{
						variantContextWriter.add(ctx2);
						this.lastVariantCount++;
						}
					}
				// check STOP File
				final long now = System.currentTimeMillis();
				if((now - lastTick) > 10*1000)//10sec
					{
					lastTick = now;
					if(STOP_FILE.exists())
						{
						LOG.warn("STOP FILE detected "+STOP_FILE+" Aborting.");
						fail_flag=true;
						break;
						}
					}
				}
			progress.finish();
			iter.close();
			iter=null;
			
			variantContextWriter.close();
			variantContextWriter=null;
			
			return outVcfFile.getPath();
			}
		catch(final Exception err)
			{
			fail_flag=true;
			LOG.error(err);
			throw new IOException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(variantContextWriter);
			if(fail_flag)
				{
				if(outVcfFile!=null && outVcfFile.exists())
					{
					LOG.warn("deleting "+outVcfFile);
					outVcfFile.delete();
					}
				if(outVcfIndexFile!=null && outVcfIndexFile.exists())
					{
					LOG.warn("deleting "+outVcfIndexFile);
					outVcfIndexFile.delete();
					}
				}
	
			}
		}
	
/** number of variant echoed the last time we called processVcfMulti */
public int getLastVariantCount() {
	return this.lastVariantCount;
	}

public static void main(String[] args) {
	System.err.println("This is a library. It's expected to be run into knime.");
	System.exit(-1);
	}
	
}
