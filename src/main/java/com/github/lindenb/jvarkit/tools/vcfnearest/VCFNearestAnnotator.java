/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfnearest;

import java.io.Closeable;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineCodec;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;
import com.github.lindenb.jvarkit.variant.infotable.VCFInfoTableModelFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VCFNearestAnnotator implements VariantAnnotator {
	private static final Logger LOG = Logger.of(VCFNearestAnnotator.class);
	private String tag = "NEARSET";
	private int max_distance = 0;
	private int limit_count = -1;
	private Path bedPath;
	private final List<String> bedColumns = new ArrayList<>();
	private Scanner scanner=null;
	private boolean use_tabix_index =false;
	
	private static abstract class Scanner implements Closeable {
		abstract void loadBed(final VCFHeader header,Path bed)throws IOException;
		abstract Stream<BedLine> query(Locatable loc)throws IOException;
		}
	private static class TreeMapScanner extends Scanner {
		private final IntervalTreeMap<List<BedLine>> loc2beds = new IntervalTreeMap<>();

		private void addInterval(final BedLine rec) {
			final Interval r=new Interval(rec);
			List<BedLine> genes = loc2beds.get(r);
			if(genes==null) {
				genes=new ArrayList<>();
				loc2beds.put(r, genes);
				}
			genes.add(rec);
			}
		
		@Override
		void loadBed(final VCFHeader header,Path bedPath) throws IOException
			{
			final SAMSequenceDictionary dict= header.getSequenceDictionary();
			final UnaryOperator<String> converter = dict==null?C->C: ContigNameConverter.fromOneDictionary(dict);
			try(BedLineReader br = new BedLineReader(bedPath)) {
				br.setContigNameConverter(converter);
				br.stream().forEach(this::addInterval);
				}
			}
		final Stream<BedLine> query(Locatable loc) throws IOException{
			return this.loc2beds.getOverlapping(loc).
				stream().
				flatMap(L->L.stream())
				;
				
			}
		@Override
		public void close() throws IOException {
			loc2beds.clear();
			}
		}
	
	private static class TabixScanner extends Scanner {
		TabixReader reader=null;
		final BedLineCodec codec=new BedLineCodec();
		UnaryOperator<String> to_tabix_ctg;
		@Override
		void loadBed(final VCFHeader header,Path bed) throws IOException
			{
			this.reader = new TabixReader(bed.toAbsolutePath().toString());
			final SAMSequenceDictionary dict= header.getSequenceDictionary();
			final UnaryOperator<String> converter = dict==null?C->C: ContigNameConverter.fromOneDictionary(dict);
			codec.setContigNameConverter(converter);
			this.to_tabix_ctg = ContigNameConverter.fromContigSet(this.reader.getChromosomes());
			}
		
		final Stream<BedLine> query(Locatable loc) throws IOException {
			final String ctg2 = this.to_tabix_ctg.apply(loc.getContig());
			if(StringUtils.isBlank(ctg2)) return Stream.empty();
			final  TabixReader.Iterator iter=reader.query(ctg2, Math.min(0,loc.getStart()-1), loc.getEnd());
			final List<BedLine>  L = new ArrayList<>();
			for(;;) {
				String line=iter.next();
				if(line==null) break;
				final BedLine rec= this.codec.decode(line);
				if(rec==null || !rec.overlaps(loc)) continue;
				L.add(rec);
				}
			return L.stream();
			}
		@Override
		public void close() throws IOException {
			if(reader!=null) try {reader.close();} catch(Throwable err) {}
			reader=null;
			}
	}
	
	public VCFNearestAnnotator() {
		}


	public void setUseTabixIndex(boolean use_tabix_index) {
		this.use_tabix_index = use_tabix_index;
		}
	public void setBedPath(Path bedPath) {
		this.bedPath = bedPath;
		}
	
	public void setColumns(final List<String> bedcolumns) {
		this.bedColumns.clear();
		this.bedColumns.addAll(bedcolumns);
	}
	
	
	@Override
	public void fillHeader(final VCFHeader header) {
		final List<String> cols2 = new ArrayList<>();
		this.bedColumns.stream().
			filter(S->!StringUtils.isBlank(S)).
			forEach(S->cols2.add(S));
		cols2.add("side");
		cols2.add("distance");
		
		
		
		header.addMetaDataLine(new VCFInfoHeaderLine(
				getTag(), 
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,"Features located less than "+this.max_distance+"bp from variant. "+
						VCFInfoTableModelFactory.encode(cols2)
				));
		
		// load BED
		if(this.bedPath!=null) {
			if(use_tabix_index && bedPath.getFileName().toString().endsWith(".bed.gz") &&
					Files.exists(Paths.get(ParsingUtils.appendToPath(bedPath.toAbsolutePath().toString(), FileExtensions.TABIX_INDEX))))
				{
				this.scanner=new TabixScanner();
				}
			else
				{
				this.scanner=new TreeMapScanner();
				}
			
			try {
				this.scanner.loadBed(header, bedPath);
				}
			catch(IOException err) {
				LOG.error(err);
				throw new RuntimeIOException(err);
				}
			}
		else
			{
			throw new IllegalArgumentException("undefined bed");
			}
		}
	
	private int distanceTo(final VariantContext ctx,final BedLine reg) {
		if(ctx.overlaps(reg)) return 0;
		else if(ctx.getEnd() < reg.getStart()) {
			return CoordMath.getLength(ctx.getEnd(), reg.getStart());
			}
		else
			{
			return CoordMath.getLength(reg.getEnd(), ctx.getStart());
			}
		}
	
	public void setTag(String tag) {
		this.tag = tag;
		}
	
	public String getTag() {
		return tag;
		}
	
	public void setMaxDistance(int max_distance) {
		this.max_distance = max_distance;
		}
	
	public void setLimitCount(int limit_count) {
		this.limit_count = limit_count;
		}
	
	@Override
	public List<VariantContext> annotate(final VariantContext ctx) throws IOException {
		final Interval r = new Interval(
				ctx.getContig(), 
				Math.max(1, ctx.getStart()-max_distance),
				ctx.getEnd()+max_distance
				);
		final List<BedLine> L1= this.scanner.query(r).
				sorted((A,B)->Integer.compare(distanceTo(ctx,A), distanceTo(ctx,B))).
				collect(Collectors.toList());
		
		if(L1.isEmpty()) return Collections.singletonList(ctx);
		final Set<String> set=new LinkedHashSet<>();
		for(BedLine rec:L1) {
			if(this.limit_count>=0 && set.size() > this.limit_count) break;
			final int side;
			if(ctx.overlaps(rec))  {
				side= 0;
				}
			else if(rec.getEnd() < ctx.getStart() ) {
				side=-1;
				}
			else
				{
				side=1;
				}
			final List<String> row = new ArrayList<>(this.bedColumns.size()+2);
			for(int x=0; x< this.bedColumns.size();++x) {
				if(StringUtils.isBlank(this.bedColumns.get(x))) continue;
				row.add(rec.getOrDefault(x,""));
				}
			row.add(String.valueOf(side));
			row.add(String.valueOf(distanceTo(ctx,rec)));
			
			set.add(String.join("|", row));
			}
		return Collections.singletonList(
				new VariantContextBuilder(ctx).
					attribute(getTag(), new ArrayList<String>(set)).
					make()
				);
		}
	
	@Override
	public void close() {
		if(scanner!=null) try {
			scanner.close();
			} catch(IOException err) {
			}
		scanner=null;
		}
}
