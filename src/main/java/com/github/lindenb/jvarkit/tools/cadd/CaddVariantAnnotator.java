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
package com.github.lindenb.jvarkit.tools.cadd;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.bio.AcidNucleics;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.ExtendedLocatable;
import com.github.lindenb.jvarkit.tabix.AbstractTabixVariantAnnotator;


public class CaddVariantAnnotator extends AbstractTabixVariantAnnotator {
	private static final Logger LOG = Logger.of(CaddVariantAnnotator.class);

	/** global can be used by vcf2r for Matilde */
	public static final String DEFAULT_CADD_FLAG_SCORE = "CADD_SCORE";
	/** global can be used by vcf2r for Matilde */
	public static final String DEFAULT_CADD_FLAG_PHRED = "CADD_PHRED";
	
	private int buffer_distance= 1_000;
	private float NA_value = -999f;
	private FileHeader fileHeader = null;
	private int column_index_for_Alt = -1;
	private int column_index_for_RawScore = -1;
	private int column_index_for_PHRED = -1;
	private Locatable lastInterval = null;
	private final List<Record> buffer = new ArrayList<>();

	
	public void setBufferDistance(int buffer_distance) {
		this.buffer_distance = Math.max(buffer_distance,10);
		}
	public void setNAValue(float NA_value) {
		this.NA_value = NA_value;
		}
	
	private class Record implements ExtendedLocatable {
		final String contig;
		final int pos;
		final Allele ref;
		final Allele alt;
		final float score;
		final float phred;
		Record(final String contig,final List<String> tokens) {
			if(tokens.size()<6) throw new JvarkitException.TokenErrors("Bad CADD line . Expected at least 6 fields",tokens);
			this.contig = contig;
			this.pos= Integer.parseInt(tokens.get(1));
			this.ref = Allele.create(tokens.get(2),true);
			this.alt = Allele.create(tokens.get(CaddVariantAnnotator.this.column_index_for_Alt),false);
			this.score = Float.parseFloat(tokens.get(CaddVariantAnnotator.this.column_index_for_RawScore));
			this.phred = Float.parseFloat(tokens.get(CaddVariantAnnotator.this.column_index_for_PHRED));
			}
		@Override
		public String getContig()
			{
			return this.contig;
			}
		@Override
		public int getStart()
			{
			return this.pos;
			}
		@Override
		public int getEnd()
			{
			return this.pos;
			}
		}
	
	public CaddVariantAnnotator(final String uri) throws IOException {
		super(uri);
		for(;;)
			{
			final String head = super.tabixReader.readLine();
			if(head==null) throw new IOException("cannot read first line of "+uri);
			if(!head.startsWith("#")) break;
			if(head.startsWith("#Chrom\tPos\tRef") ||
			   head.startsWith("#Chr\tPos\tRef"))
				{
				this.fileHeader = new FileHeader(head, S->Arrays.asList(CharSplitter.TAB.split(S)));
				break;
				}
			}
		if(this.fileHeader==null || this.fileHeader.isEmpty())
			{
			throw new IOException("Cannot read tabix file header");
			}
		
		if(!(this.fileHeader.get(0).equals("#Chrom") || this.fileHeader.get(0).equals("#Chr"))) {
			throw new IOException("Illegal Header: Cannot get column index of #Chrom at 0");
			}
		this.fileHeader.assertColumn("Pos",1);
		this.fileHeader.assertColumn("Ref",2);
		this.column_index_for_Alt = this.fileHeader.getColumnIndex("Alt");
		this.column_index_for_RawScore = this.fileHeader.getColumnIndex("RawScore");
		this.column_index_for_PHRED = this.fileHeader.getColumnIndex("PHRED");
		}
	
	
	private List<Record> fetch(final VariantContext loc) {
		final String chromCadd = this.contig(loc);
		if(StringUtils.isBlank(chromCadd)) return Collections.emptyList();
		final Locatable query = new SimpleInterval(chromCadd,loc.getStart(),loc.getEnd());
		if(this.lastInterval==null || !this.lastInterval.contains(query)) {
			this.buffer.clear();
			this.lastInterval = new SimpleInterval(chromCadd, Math.max(0, loc.getStart()-1), Math.max(loc.getEnd(), loc.getStart()+this.buffer_distance));
				try {
					final TabixReader.Iterator iter = super.tabixReader.query(
				
						this.lastInterval.getContig(),
						this.lastInterval.getStart(),
						this.lastInterval.getEnd()
						);
				for(;;) {
					final String line=iter.next();
					if(line==null) break;
					if(line.startsWith("#") || StringUtil.isBlank(line))continue;
					final List<String> tokens = this.fileHeader.split(line);
				
					if(!AcidNucleics.isATGC(tokens.get(2)))
						{
						LOG.warn("REF allele not suitable in  line "+line+". skipping");
						continue;
						}
					if(!AcidNucleics.isATGC(tokens.get(this.column_index_for_Alt)))
						{
						LOG.warn("ALT allele not suitable in  line "+line+". skipping");
						continue;
						}
				
					if(!tokens.get(0).equals(chromCadd)) throw new IllegalStateException("Expected CADD contig "+chromCadd+" in "+line);
	
					final Record rec=new Record(loc.getContig(),tokens);
	
					this.buffer.add(rec);
					}
				} 
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		return this.buffer.stream().filter(R->
				R.getStart() == loc.getStart() &&
				R.ref.equals(loc.getReference()) &&
				loc.getAlternateAlleles().contains(R.alt)
				).collect(Collectors.toList());
		}

	private void  checkRange(float f) {
		if(f <= this.NA_value) {
			throw new IllegalArgumentException("Got score = "+f +" but the value for N/A is "+this.NA_value);
			}
		}
	
	@Override
	public List<VariantContext> annotate(final VariantContext ctx) throws IOException {
			return Collections.singletonList(annotate1(ctx));
			}
	
	private VariantContext annotate1(final VariantContext ctx) throws IOException {
			final List<Record> cadd_rec_for_ctx = fetch(ctx);
			if(cadd_rec_for_ctx.isEmpty()) return ctx;
			
			final List<Float> cadd_array_score=new ArrayList<>();
			final List<Float> cadd_array_phred=new ArrayList<>();
			boolean got_non_null = false;
			for(final Allele alt:ctx.getAlternateAlleles()) {
				final Record rec = cadd_rec_for_ctx.
						stream().
						filter(REC->REC.alt.equals(alt)).
						findAny().
						orElse(null);
				if(rec==null) {
					cadd_array_score.add(NA_value);
					cadd_array_phred.add(NA_value);
					}
				else {
					got_non_null = true;
					checkRange(rec.score);
					checkRange(rec.phred);
					cadd_array_score.add(rec.score);
					cadd_array_phred.add(rec.phred);
					}
				}

			if(!cadd_array_score.isEmpty() && got_non_null)
				{
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.attribute(DEFAULT_CADD_FLAG_SCORE, cadd_array_score);
				vcb.attribute(DEFAULT_CADD_FLAG_PHRED, cadd_array_phred);
				return vcb.make();
				}
		return ctx;
		}
	

	@Override
	public void fillHeader(final VCFHeader header) {
		
		if(header.getInfoHeaderLine(DEFAULT_CADD_FLAG_PHRED)!=null) {
			throw new JvarkitException.DuplicateVcfHeaderInfo(header, DEFAULT_CADD_FLAG_PHRED);
		}
		if(header.getInfoHeaderLine(DEFAULT_CADD_FLAG_SCORE)!=null) {
			throw new JvarkitException.DuplicateVcfHeaderInfo(header, DEFAULT_CADD_FLAG_SCORE);
		}

		
		header.addMetaDataLine(new VCFInfoHeaderLine(
			DEFAULT_CADD_FLAG_SCORE,
			VCFHeaderLineCount.A,
			VCFHeaderLineType.Float,
			"Score suggests that that variant is likely to be  observed (negative values) vs simulated(positive values)."+
			"However, raw values do have relative meaning, with higher values indicating that a variant is more likely to be simulated (or -not observed-) and therefore more likely to have deleterious effects." +
			" URI was " +getUri() + 
			". We use "+NA_value+" for unknown value."
			));
		header.addMetaDataLine(new VCFInfoHeaderLine(
			DEFAULT_CADD_FLAG_PHRED,
			VCFHeaderLineCount.A,
			VCFHeaderLineType.Float,
			"PHRED expressing the rank in order of magnitude terms. For example, reference genome single nucleotide variants at the 10th-% of CADD scores are assigned to CADD-10, top 1% to CADD-20, top 0.1% to CADD-30, etc. " +
			" URI was " +getUri() + 
			". We use "+NA_value+" for unknown value."
			));
		}
	
	}
