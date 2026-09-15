/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcfalphagenome;

import java.io.Closeable;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
 
 BEGIN_DOC
 
 ## Example
 
 ```
 $ bcftools view in.vcf |\
   java -jar dist/jvarkit.jar vcfalphagenome \
 	--impact   /path/to/google/AlphaGenome/alphagenome_variant_impact_score_snvs.GRCh38.tsv.gz \\
 	--splicing /path/to/google/AlphaGenome/combined_alphagenome_splicing_snvs.GRCh38.tsv.gz |\
 	grep -v "##" | head -n 4

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
chr1    26169987        .       G       A       .       .       AG_PHRED=22.68624;AG_RAW=0.7722;AG_SPLICING=0.08169
chr1    42759328        .       G       A       .       .       AG_PHRED=18.73157;AG_RAW=0.4941;AG_SPLICING=0.06523
chr1    43558671        .       C       T       .       .       AG_PHRED=14.63943;AG_RAW=0.307;AG_SPLICING=0.05377
```
 
 END_DOC
 
 */

@Program(
		name="vcfalphagenome",
		description= "Annotate VCF with Tabix index files from Google AlphaGenome",
		creationDate="20260915",
		modificationDate="20260915",
		keywords={"vcf","prediction","google","alphagenome"},
		jvarkit_amalgamion = true,
		menu="VCF Manipulation"
		)
public class VcfAlphaGenome extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfAlphaGenome.class);
	@Parameter(names={"-d","--buffer-size"},description="Buffer size / processing window size. " + DistanceParser.OPT_DESCRIPTION , converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int buffer_distance= 1_000;
	@Parameter(names={"--impact"},description="Path to tabix indexed impact_score_snvs.tsv.gz")
	private String impactPath=null;
	@Parameter(names={"--splicing","--combined"},description="Path to tabix indexed combined_splicing_snvs.tsv.gz")
	private String splicingPath=null;
	@Parameter(names={"-na","--na"},description="value  used for 'allele-not-found'.")
	private float NA_value = -999f;
	@Parameter(names={"--tag"},description="INFO tag prefix.")
	private String tag_prefix = "AG_";

	private static abstract class Record {
		@SuppressWarnings("unused")
		String contig;
		final int pos;
		final Allele ref;
		final Allele alt;
		Record(final String[] tokens) {
			this.contig = tokens[0];
			this.pos = Integer.parseInt(tokens[1]);
			this.ref = Allele.create(tokens[2], true);	
			this.alt = Allele.create(tokens[3], false);	
			}
		}
	
	private static class ImpactRecord extends Record {
		final float raw_score;
		final float phred;
		
		ImpactRecord(final String[] tokens) {
			super(tokens);
			this.raw_score = Float.parseFloat(tokens[4]);
			this.phred = Float.parseFloat(tokens[5]);
			}
		
		}
	private static class SplicingRecord extends Record {
		final float alphagenome_splicing;
		SplicingRecord(final String[] tokens) {
			super(tokens);
			this.alphagenome_splicing = Float.parseFloat(tokens[4]);
			}
	}
	
	private  abstract class DataSource<REC extends Record> implements Closeable {
		final TabixReader tabix;
		private Locatable lastInterval = null;
		private final List<REC> buffer = new ArrayList<>();
		private final ContigNameConverter converter;

		protected DataSource(final String uri,final String expect_header) throws IOException {
			this.tabix = new TabixReader(uri);
			this.converter = ContigNameConverter.fromContigSet(tabix.getChromosomes());
			final String head = this.tabix.readLine();
			if(head==null) throw new IOException("cannot read first line of "+uri);
			if(!head.startsWith(expect_header)) {
				throw new IOException("in "+uri+" expected header: "+expect_header+" but got "+head);
				}
							
			}
		
		protected void  checkRange(final float f) {
			if(f <= VcfAlphaGenome.this.NA_value) {
				throw new IllegalArgumentException("Got score = "+f +" but the value for N/A is "+ VcfAlphaGenome.this.NA_value);
				}
			}

		
		/** convert variant contig */
		protected String contig(final VariantContext ctx) {
			return converter.apply(ctx.getContig());
			}

		abstract void fillHeader(final VCFHeader header);
		abstract void annotate(final VariantContext ctx,final VariantContextBuilder vcb);
		protected abstract REC makeRecord(String[] tokens);
		
		protected List<REC> fetch(final VariantContext loc) {
			final String chromAlphaGenome = this.contig(loc);
			if(StringUtils.isBlank(chromAlphaGenome)) return Collections.emptyList();
			final Locatable query = new SimpleInterval(chromAlphaGenome,loc.getStart(),loc.getEnd());
			if(this.lastInterval==null || !this.lastInterval.contains(query)) {
				this.buffer.clear();
				this.lastInterval = new SimpleInterval(chromAlphaGenome, Math.max(0, loc.getStart()-1), Math.max(loc.getEnd(), loc.getStart()+ VcfAlphaGenome.this.buffer_distance));
					try {
						final TabixReader.Iterator iter = this.tabix.query(
					
							this.lastInterval.getContig(),
							this.lastInterval.getStart(),
							this.lastInterval.getEnd()
							);
					for(;;) {
						final String line=iter.next();
						if(line==null) break;
						if(line.startsWith("#") || StringUtil.isBlank(line))continue;
						final String[] tokens = CharSplitter.TAB.split(line);
						final REC rec = makeRecord(tokens);
		
						rec.contig = loc.getContig();
		
						this.buffer.add(rec);
						}
					} 
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
				}
			return this.buffer.stream().filter(R->
					R.pos == loc.getStart() &&
					R.ref.equals(loc.getReference()) &&
					loc.getAlternateAlleles().contains(R.alt)
					).collect(Collectors.toList());
			}

		
		@Override
		public void close()  {
			this.tabix.close();
			}
		}
	private class ImpactDataSource extends DataSource<ImpactRecord> {
		final VCFInfoHeaderLine info_raw;
		final VCFInfoHeaderLine info_phred;
		ImpactDataSource(final String uri)  throws IOException {
			super(uri,"#CHROM\tPOS\tREF\tALT\traw_score\tPHRED");
			this.info_raw = new VCFInfoHeaderLine(
					VcfAlphaGenome.this.tag_prefix+ "RAW",
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"Alpha genome raw score. URI was : "+uri + "We use "+NA_value+" for unknown value."
					);
			this.info_phred = new VCFInfoHeaderLine(
					VcfAlphaGenome.this.tag_prefix+ "PHRED",
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"Alpha genome PHRED score. URI was : "+uri + "We use "+NA_value+" for unknown value. Phred scores range from 0 (no impact) to 70 (maximal impact)"
					);
			}
		@Override
		protected ImpactRecord makeRecord(final String[] tokens) {
			return new ImpactRecord(tokens);
			}
		@Override
		public void fillHeader(final VCFHeader header) {
			header.addMetaDataLine(this.info_raw);
			header.addMetaDataLine(this.info_phred);
			}
		@Override
		void annotate(VariantContext ctx,final VariantContextBuilder vcb) {
			final List<ImpactRecord> rec_for_ctx = fetch(ctx);
			if(rec_for_ctx.isEmpty()) return;
			
			final List<Float> array_score=new ArrayList<>();
			final List<Float> array_phred=new ArrayList<>();
			boolean got_non_null = false;
			for(final Allele alt:ctx.getAlternateAlleles()) {
				final ImpactRecord rec = rec_for_ctx.
						stream().
						filter(REC->REC.alt.equals(alt)).
						findAny().
						orElse(null);
				if(rec==null) {
					array_score.add(NA_value);
					array_phred.add(NA_value);
					}
				else {
					got_non_null = true;
					checkRange(rec.raw_score);
					checkRange(rec.phred);
					array_score.add(rec.raw_score);
					array_phred.add(rec.phred);
					}
				}
	
			if(!array_phred.isEmpty() && got_non_null)
				{
				vcb.attribute(this.info_raw.getID(), array_score);
				vcb.attribute(this.info_phred.getID(), array_phred);
				}
			}
		}

	private class CombinedSplicingDataSource extends DataSource<SplicingRecord> {
		final VCFInfoHeaderLine info_splicing;
		CombinedSplicingDataSource(final String uri)  throws IOException {
			super(uri,"#CHROM\tPOS\tREF\tALT\talphagenome_splicing");
			this.info_splicing = new VCFInfoHeaderLine(
					VcfAlphaGenome.this.tag_prefix + "SPLICING",
					VCFHeaderLineCount.A,
					VCFHeaderLineType.Float,
					"Alpha genome splicing score. URI was : "+uri + "We use "+NA_value+" for unknown value. score=0 meaning no impact on splicing and score>0 meaning variant is predicted to alter splice site usage. A higher value indicates a stronger disruptive effect on splice site usage (a splicing-disruptive effect)"
					);
			}
		@Override
		protected SplicingRecord makeRecord(final String[] tokens) {
			return new SplicingRecord(tokens);
			}
		@Override
		public void fillHeader(final VCFHeader header) {
			header.addMetaDataLine(this.info_splicing);
			}
		@Override
		void annotate(final VariantContext ctx,final VariantContextBuilder vcb) {
			final List<SplicingRecord> rec_for_ctx = fetch(ctx);
			if(rec_for_ctx.isEmpty()) return;
			
			final List<Float> array_score=new ArrayList<>();
			boolean got_non_null = false;
			for(final Allele alt:ctx.getAlternateAlleles()) {
				final SplicingRecord rec = rec_for_ctx.
						stream().
						filter(REC->REC.alt.equals(alt)).
						findAny().
						orElse(null);
				if(rec==null) {
					array_score.add(NA_value);
					}
				else {
					got_non_null = true;
					checkRange(rec.alphagenome_splicing);
					array_score.add(rec.alphagenome_splicing);
					}
				}
	
			if(!array_score.isEmpty() && got_non_null)
				{
				vcb.attribute(this.info_splicing.getID(), array_score);
				}
			}

		}

	
	private final List<DataSource<?>> dataSources = new ArrayList<>();
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeVcf() {
		
		try {
			if(!StringUtils.isBlank(this.impactPath)) {
				this.dataSources.add(new ImpactDataSource(this.impactPath));
				}
			if(!StringUtils.isBlank(this.splicingPath)) {
				this.dataSources.add(new CombinedSplicingDataSource(this.splicingPath));
				}
			}
		catch(final Throwable err) {
			getLogger().error(err);
			return -1;
			}
		if(this.dataSources.isEmpty()) {
			getLogger().error("both input impact/splicing are blank");
			return -1;
			}
		return super.beforeVcf();
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		final VCFHeader h2 = new VCFHeader(iterin.getHeader());
		for(DataSource<?> ds: this.dataSources) {
			ds.fillHeader(h2);
			}
		JVarkitVersion.getInstance().addMetaData(this, h2);
		out.writeHeader(h2);
		while(iterin.hasNext()) {
			final VariantContext ctx=iterin.next();
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			for(DataSource<?> ds: this.dataSources) {
				ds.annotate(ctx,vcb);
				}
			out.add(vcb.make());
			}
		
		return 0;
		}
	
	@Override
	protected void afterVcf() {
		for(DataSource<?> ds: this.dataSources) {
			ds.close();
			}
		super.afterVcf();
		}

	
	public static void main(final String[] args) {
		new VcfAlphaGenome().instanceMainWithExit(args);

	}

}
