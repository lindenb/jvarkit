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
package com.github.lindenb.jvarkit.tools.evai;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

For Adeline Goudal (master Bioinfo 2019) .

END_DOC
*/

@Program(name="vcfevaiannot",
	description="Annotate vcf with evai data.",
	keywords= {"vcf","annotation","evai"},
	creationDate = "20190507",
	modificationDate="20190507",
	generate_doc = false
	)
public class VcfEvaiAnnot extends Launcher {


	private static final Logger LOG = Logger.build(VcfEvaiAnnot.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-T","--tabix"},description="File containing path to eVAI files. One file per line. Each file MUST be compressed and indexed with tabix.",required = true)
	private Path tabixListPath = null;
	@Parameter(names={"-B","--buffer-size"},description="Buffer size. "+DistanceParser.OPT_DESCRIPTION ,converter = DistanceParser.StringConverter.class,splitter = NoSplitter.class)
	private int buffer_size = 1_000;
	@Parameter(names={"--prefix"},description="Attribute prefix")
	private String prefix = "EVAI_";
	@Parameter(names={"--ignore-filtered"},description="Ignore FILTERed variants (should go faster)")
	private boolean ignore_filtered = false;

	
	private class EvaiLine implements Locatable {
		final String tokens[];
		EvaiLine(final String tokens[]) {
			this.tokens = tokens;
			}
		@Override
		public String getContig() {
			return tokens[0];
			}
		
		@Override
		public int getStart() {
			return Integer.parseInt(tokens[1]);
			}
		@Override
		public int getEnd() {
			return Integer.parseInt(tokens[2]);
			}
		
		@Override
		public String toString() {
			return String.join("\t", tokens);
			}
		
		boolean match(final VariantContext ctx) {
			if(!this.getContig().equals(ctx.getContig())) return false;
			if(this.getStart()!=ctx.getStart()) return false;
			if(this.getEnd()!=ctx.getEnd()) return false;
			if(!ctx.getReference().getDisplayString().equals(tokens[3])) return false;
			if(ctx.getNAlleles()==1) return false;
			if(!ctx.getAlleles().get(1).getDisplayString().equals(tokens[4])) return false;
			return true;
			}
		}
	
	private class EvaiTabix implements Closeable {
		final String filename;
		final TabixReader tabix;
		String sample = null;
		Interval lastInterval = null;
		final List<EvaiLine> buffer = new ArrayList<>();
		final Map<String,Integer> column2index = new HashMap<>();
		
		EvaiTabix(final String filename) throws IOException {
			this.filename = filename;
			this.tabix = new TabixReader(filename);
			final String version_tag = "##eVAI-version=";//0.4.2
			final String sample_tag = "##SAMPLE-ID:";
			for(;;) {
				String line = this.tabix.readLine();
				if(line==null || !line.startsWith("#")) break;
				if(line.startsWith(version_tag)) {
					final String version = line.substring(version_tag.length()).trim();
					if(!version.equals("0.4.2")) throw new IOException("unknown vai version "+ line+" "+filename);
					}
				if(line.startsWith(version_tag)) {
					final String version = line.substring(version_tag.length()).trim();
					if(!version.equals("0.4.2")) throw new IOException("unknown vai version "+ line+" "+filename);
					}
				if(line.startsWith(sample_tag)) {
					this.sample = line.substring(sample_tag.length()).trim();
					}
				if(line.startsWith("#CHROM")) {
					final String tokens[] = CharSplitter.TAB.split(line);
					for(int i=0;i< tokens.length;++i) {
						this.column2index.put(tokens[i],i);
						}
					}
				}
			if(StringUtils.isBlank(this.sample)) {
				throw new IOException("No sample found in "+this.filename);
				}
			if(this.column2index.isEmpty()) {
				throw new IOException("No header found in "+this.filename);
				}
			}
		
		
		
		Optional<EvaiLine> query(final VariantContext ctx) {
			if(this.lastInterval==null ||
				!this.lastInterval.getContig().equals(ctx.getContig()) ||
				!CoordMath.encloses(lastInterval.getStart(), lastInterval.getEnd(),
						ctx.getStart(),ctx.getEnd())) {
				
				buffer.clear();
				this.lastInterval = new Interval(ctx.getContig(),
						Math.max(1,ctx.getStart()-10),
						ctx.getEnd() + Math.max(1, buffer_size)
						);
				final TabixReader.Iterator iter = this.tabix.query(this.lastInterval.getContig(), this.lastInterval.getStart(), this.lastInterval.getEnd());
				final CharSplitter tab = CharSplitter.TAB;
				try {
					for(;;) {
						final String line = iter.next();
						if(line==null) break;
						final EvaiLine evai = new EvaiLine(tab.split(line));
						this.buffer.add(evai);
						}
					}
				catch(final IOException err) {
					throw new RuntimeIOException(err);
					}
				}
			return this.buffer.stream().filter(L->L.match(ctx)).findFirst();
			}
		
		@Override
		public void close(){
			this.tabix.close();
			}
		}

	private final Map<String,EvaiTabix> sample2tabix = new HashMap<>();
	
	private boolean isBooleanField(final String T) {
		return T.startsWith("BP") || T.startsWith("BS") || T.startsWith("PM") || T.startsWith("PP") || T.startsWith("PS")|| T.startsWith("PV");
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		final VCFHeader header0 = iterin.getHeader();
		final VCFHeader header2 = new VCFHeader(header0);
		
		final Set<VCFFormatHeaderLine> meta = new HashSet<>();
		this.sample2tabix.values().
			stream().
			flatMap(T->T.column2index.keySet().stream()).
			filter(T->isBooleanField(T)).
			map(T->new VCFFormatHeaderLine(this.prefix+T,1,VCFHeaderLineType.Integer,T)).
			forEach(H->meta.add(H));
		
		meta.add(new VCFFormatHeaderLine(this.prefix+"FINAL_CLASSIFICATION",1,VCFHeaderLineType.String,"FINAL_CLASSIFICATION"));
		
		meta.stream().forEach(M->header2.addMetaDataLine(M));
		
		JVarkitVersion.getInstance().addMetaData(this, header2);
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header0).logger(LOG).build();
		
		out.writeHeader(header2);
		while(iterin.hasNext()) {
			final VariantContext ctx = progress.apply(iterin.next());
			if(this.ignore_filtered && ctx.isFiltered()) {
				out.add(ctx);
				continue;
				}
			
			final List<Genotype> genotypes = new ArrayList<>(ctx.getNSamples());
			for(final Genotype g: ctx.getGenotypes()) {
				if(g.isHomRef() || g.isNoCall()) {
					genotypes.add(g);
					continue;
					}

				final EvaiTabix tbx = this.sample2tabix.get(g.getSampleName());
				if(tbx==null) {
					genotypes.add(g);
					continue;
					}
				final Optional<EvaiLine> candidate = tbx.query(ctx);
				if(!candidate.isPresent()) {
					genotypes.add(g);
					continue;
					}
				final EvaiLine line = candidate.get();
				final GenotypeBuilder gb= new GenotypeBuilder(g);
				for(final String column : tbx.column2index.keySet()) {
					if(!isBooleanField(column)) continue;
					final int col_idx = tbx.column2index.get(column);
					if(col_idx>= line.tokens.length) continue;
					final String valuestr=line.tokens[col_idx];
					if(valuestr.equals("n.a.")) continue;
					final Object value ;
					if(valuestr.equals("true")) value=1;
					else if(valuestr.equals("TRUE(STRONG)")) value=10;
					else if(valuestr.equals("false")) value=0;
					else throw new IllegalArgumentException(column+" "+valuestr);
					gb.attribute(this.prefix+column, value);
					}
				genotypes.add(gb.make());
				}
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.genotypes(genotypes);
			out.add(vcb.make());
			}
		
		progress.close();
		return 0;
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			Files.lines(this.tabixListPath).forEach(L->{
				try {
					final EvaiTabix t = new EvaiTabix(L);
					this.sample2tabix.put(t.sample, t);
					}
				catch(final IOException err) {
					throw new RuntimeIOException(err);
					}
				});
			
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			this.sample2tabix.values().stream().forEach(T->T.close());
			}
		}
	
	public static void main(final String[] args) {
		new VcfEvaiAnnot().instanceMainWithExit(args);
		}
}
