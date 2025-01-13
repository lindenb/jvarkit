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
package com.github.lindenb.jvarkit.tools.vcfpolyx;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.ChromosomeSequence;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.variant.VariantAnnotator;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class PolyXVariantAnnotator implements VariantAnnotator {
	private OptionalInt filterTreshold = OptionalInt.empty();
	private String polyXtag = "POLYX";
	private boolean skipFiltered=false;
	private final ReferenceSequenceFile  referenceSequenceFile;
	private VCFInfoHeaderLine infoHeaderLine = null;
	private VCFFilterHeaderLine filterHeaderLine = null;
	private ChromosomeSequence genomicContig=null;
	private final ContigNameConverter contigNameConverter;

	
	public PolyXVariantAnnotator(Path faixPath) throws IOException {
		this.referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(faixPath);
		this.contigNameConverter = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(referenceSequenceFile));
		}

	public void setSkipFiltered(boolean skipFiltered) {
		this.skipFiltered = skipFiltered;
	}
	
	public boolean isSkipFiltered() {
		return skipFiltered;
	}
	public void setTag(String polyXtag) {
		this.polyXtag = polyXtag;
		if(StringUtils.isBlank(polyXtag)) throw new IllegalArgumentException("tag cannot be empty");
	}
	public String getTag() {
		return polyXtag;
	}
	
	public void setMaxRepeat(OptionalInt filterTreshold) {
		this.filterTreshold = filterTreshold;
		}
	
	public OptionalInt getMaxRepeat() {
		return filterTreshold;
		}
	
	@Override
	public void fillHeader(VCFHeader header) {
		
		this.infoHeaderLine = new VCFInfoHeaderLine(
				getTag(),
				1,
				VCFHeaderLineType.Integer,
				"Number of repeated bases around REF")
				;
		header.addMetaDataLine(this.infoHeaderLine);

		
		
		if(getMaxRepeat().isPresent()) {
			this.filterHeaderLine = new VCFFilterHeaderLine(
					infoHeaderLine.getID()+"_ge_"+ getMaxRepeat().getAsInt(),
					"Number of repeated bases around REF is greater or equal to " + getMaxRepeat().getAsInt()
					);
			header.addMetaDataLine(filterHeaderLine);
			}
		else
			{
			this.filterHeaderLine = null;
			}
		}
	
	
	@Override
	public List<VariantContext> annotate(VariantContext ctx) throws IOException {
			if(isSkipFiltered() && ctx.isFiltered())
				{
				return Collections.singletonList(ctx);
				}
		
			final String normalizedContig = this.contigNameConverter.apply(ctx.getContig());
			if(StringUtils.isBlank(normalizedContig)) {
				return Collections.singletonList(ctx);
				}
		
			if(this.genomicContig==null || !this.genomicContig.hasName(normalizedContig))
				{
				this.genomicContig= new GenomicSequence(this.referenceSequenceFile, normalizedContig);
				}
		
		final VariantContextBuilder b = new VariantContextBuilder(ctx);
	
		// https://github.com/lindenb/jvarkit/issues/165
		final boolean indel_flag = ctx.isIndel();
		
		int count=1;
		int pos0 = ctx.getStart()-1;
		// https://github.com/lindenb/jvarkit/issues/165
		if(indel_flag) {
			pos0++;
			}
		char c0 = Character.toUpperCase(genomicContig.charAt(pos0));
		//go left
		pos0--;
		while(pos0>=0 && c0==Character.toUpperCase(genomicContig.charAt(pos0)))
			{
			++count;
			pos0--;
			}
		//go right
		pos0 = ctx.getEnd()-1;
		// https://github.com/lindenb/jvarkit/issues/165
		if(indel_flag) {
			pos0++;
			}
		
		c0 = Character.toUpperCase(genomicContig.charAt(pos0));
		pos0++;
		while(pos0< genomicContig.length()
			&& c0==Character.toUpperCase(genomicContig.charAt(pos0)))
			{
			++count;
			++pos0;
			}
		b.attribute(infoHeaderLine.getID(),count);
		
		/* filter */
		if(getMaxRepeat().isPresent())
			{
			if(count >= getMaxRepeat().getAsInt()) {
				b.filter(this.filterHeaderLine.getID());
				}
			else if(!ctx.isFiltered()) {
				b.passFilters();
				}
			}
		
		return Collections.singletonList(b.make());
		}
	
	@Override
	public void close() {
		this.genomicContig = null;
		try {
			this.referenceSequenceFile.close();
			}
		catch(IOException err) {
			}
		}
}
