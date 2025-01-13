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

package com.github.lindenb.jvarkit.tabix;

import java.io.IOException;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;

/** abstract Tabix Reader to annotate Variants */
public abstract class AbstractTabixVariantAnnotator implements VariantAnnotator {
	private final String uri;
	protected final TabixReader tabixReader;
	private final ContigNameConverter converter;
	/** open tabix reader, uri can be emty or null, the annotation will be disabled */
	public AbstractTabixVariantAnnotator(final String uri) throws IOException {
		this.tabixReader = StringUtils.isBlank(uri)?null:new TabixReader(uri);
		this.uri = uri;
		this.converter = this.tabixReader==null?null:ContigNameConverter.fromContigSet(this.tabixReader.getChromosomes());
		}
	
	/** return true of variant and coordinate both overlap at 'fraction' */
	protected boolean overlaps(final VariantContext ctx,final int svStart, final int svEnd, final double fraction) {
		if(!CoordMath.overlaps(ctx.getStart(), ctx.getEnd(), svStart, svEnd)) return false;
		final int x0 = Math.max(ctx.getStart(),svStart);
		final int x1 = Math.min(ctx.getEnd(),svEnd);
		final double shared_len = CoordMath.getLength(x0, x1);
		if(shared_len/ctx.getLengthOnReference() <  fraction) return false;
		if(shared_len/CoordMath.getLength(svStart, svEnd) < fraction) return false;
		return true;
		}
	
	/** return true if tabix reader was opened */
	public boolean isValid() {
		return this.tabixReader!=null;
		}
	
	/** get filename / uri for this tabix */
	public String getUri() {
		return uri;
		}
	/** close tabix reader if any */
	@Override
	public void close() {
		if(this.tabixReader!=null) this.tabixReader.close();
		}
	/** convert variant contig */
	protected String contig(final VariantContext ctx) {
		return converter.apply(ctx.getContig());
		}
	/** return true if variant contig can be processed */
	protected boolean hasContig(final VariantContext ctx) {
		return !StringUtil.isBlank(contig(ctx));
		}
	
	@Override
	public String toString() {
		return getClass().getName()+":"+this.uri;
		}
	}
