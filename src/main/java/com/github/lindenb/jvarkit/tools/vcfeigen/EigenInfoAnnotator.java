/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
* 2016 creation

*/
package com.github.lindenb.jvarkit.tools.vcfeigen;

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
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.vcf.InfoAnnotator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class EigenInfoAnnotator implements InfoAnnotator
	{
	private static abstract class AbstractFeature implements Feature
		{
		final String contig;
		final int pos;
		final Allele ref;
		final Allele alt;
		final String tokens[];
		private AbstractFeature(final String tokens[]) {
			this.contig = tokens[0];
			this.pos = Integer.parseInt(tokens[1]);
			this.ref = Allele.create(tokens[2], true);
			this.alt = Allele.create(tokens[3], true);
			this.tokens=tokens;
			}
		
		@Override
		public String getChr()
			{
			return contig;
			}
		@Override
		public String getContig()
			{
			return contig;
			}
		@Override
		public int getStart()
			{
			return pos;
			}
		@Override
		public int getEnd()
			{
			return pos;
			}
		
		}
	private  static class CodingFeature extends AbstractFeature
		{
		private CodingFeature(final String tokens[]) {
			super(tokens);
			}
		}
	private  static class NonCodingFeature extends AbstractFeature
		{
		private NonCodingFeature(final String tokens[]) {
			super(tokens);
			}
		}
	
	
	
	private static class CodingFeatureCodec extends AsciiFeatureCodec<CodingFeature>
		{
		final Pattern tab=Pattern.compile("[\t]");
		CodingFeatureCodec() {
			super(CodingFeature.class);
			}
		@Override
		public CodingFeature decode(final String line)
			{
			if(line.startsWith("chr")) return null;
			return new CodingFeature(tab.split(line));
			}
		
		@Override
		public boolean canDecode(final String path)
			{
			return path.endsWith(".tab.gz");
			}
		@Override
		public Object readActualHeader(LineIterator arg0)
			{
			return null;
			}
		}
	private static class NonCodingFeatureCodec extends AsciiFeatureCodec<NonCodingFeature>
		{
		final Pattern tab=Pattern.compile("[\t]");
		NonCodingFeatureCodec() {
			super(NonCodingFeature.class);
			}
		@Override
		public NonCodingFeature decode(String line)
			{
			if(line.startsWith("chr")) return null;
			return new NonCodingFeature(tab.split(line));
			}
		
		@Override
		public boolean canDecode(final String path)
			{
			return path.endsWith(".tab.gz");
			}
		@Override
		public Object readActualHeader(LineIterator arg0)
			{
			return null;
			}
		}
	private File eigenDirectory = null;
	private final List<VCFInfoHeaderLine> noncodingheaderlines = new ArrayList<>();
	private TabixFeatureReader<NonCodingFeature, Object> nonCodingFeatureReader = null;
	private TabixFeatureReader<CodingFeature, Object> codingFeatureReader = null;
	private int prev_contig=-1;
		public EigenInfoAnnotator(final File dir) {
		this.eigenDirectory = dir;
		final String prefix="";
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"chr", VCFHeaderLineCount.A, VCFHeaderLineType.String, "chr"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"position", VCFHeaderLineCount.A, VCFHeaderLineType.String, "position"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"ref", VCFHeaderLineCount.A, VCFHeaderLineType.String, "ref"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"alt", VCFHeaderLineCount.A, VCFHeaderLineType.String, "alt"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"GERP_NR", VCFHeaderLineCount.A, VCFHeaderLineType.String, "GERP_NR"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"GERP_RS", VCFHeaderLineCount.A, VCFHeaderLineType.String, "GERP_RS"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"PhyloPri", VCFHeaderLineCount.A, VCFHeaderLineType.String, "PhyloPri"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"PhyloPla", VCFHeaderLineCount.A, VCFHeaderLineType.String, "PhyloPla"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"PhyloVer", VCFHeaderLineCount.A, VCFHeaderLineType.String, "PhyloVer"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"PhastPri", VCFHeaderLineCount.A, VCFHeaderLineType.String, "PhastPri"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"PhastPla", VCFHeaderLineCount.A, VCFHeaderLineType.String, "PhastPla"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"PhastVer", VCFHeaderLineCount.A, VCFHeaderLineType.String, "PhastVer"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"H3K4Me1", VCFHeaderLineCount.A, VCFHeaderLineType.String, "H3K4Me1"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"H3K4Me3", VCFHeaderLineCount.A, VCFHeaderLineType.String, "H3K4Me3"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"H3K27ac", VCFHeaderLineCount.A, VCFHeaderLineType.String, "H3K27ac"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"TFBS_max", VCFHeaderLineCount.A, VCFHeaderLineType.String, "TFBS_max"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"TFBS_sum", VCFHeaderLineCount.A, VCFHeaderLineType.String, "TFBS_sum"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"TFBS_num", VCFHeaderLineCount.A, VCFHeaderLineType.String, "TFBS_num"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"OCPval", VCFHeaderLineCount.A, VCFHeaderLineType.String, "OCPval"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"DnaseSig", VCFHeaderLineCount.A, VCFHeaderLineType.String, "DnaseSig"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"DnasePval", VCFHeaderLineCount.A, VCFHeaderLineType.String, "DnasePval"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"FaireSig", VCFHeaderLineCount.A, VCFHeaderLineType.String, "FaireSig"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"FairePval", VCFHeaderLineCount.A, VCFHeaderLineType.String, "FairePval"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"PolIISig", VCFHeaderLineCount.A, VCFHeaderLineType.String, "PolIISig"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"PolIIPval", VCFHeaderLineCount.A, VCFHeaderLineType.String, "PolIIPval"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"ctcfSig", VCFHeaderLineCount.A, VCFHeaderLineType.String, "ctcfSig"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"ctcfPval", VCFHeaderLineCount.A, VCFHeaderLineType.String, "ctcfPval"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"cmycSig", VCFHeaderLineCount.A, VCFHeaderLineType.String, "cmycSig"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"cmycPval", VCFHeaderLineCount.A, VCFHeaderLineType.String, "cmycPval"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"Eigen_raw", VCFHeaderLineCount.A, VCFHeaderLineType.String, "Eigen_raw"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"Eigen_phred", VCFHeaderLineCount.A, VCFHeaderLineType.String, "Eigen_phred"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"Eigen_PC_raw", VCFHeaderLineCount.A, VCFHeaderLineType.String, "Eigen_PC_raw"));
		this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"Eigen_PC_phred", VCFHeaderLineCount.A, VCFHeaderLineType.String, "Eigen_PC_phred"));		
			}
	
	private static int contig2eigen(String chr) {
		if(chr.toLowerCase().startsWith("chr")) chr=chr.substring(3);
		try {
			int c=Integer.parseInt(chr.trim());
			if(c<1 || c>22) return -1;
			return c;
			}
		catch(NumberFormatException err) {
			return -1;
			}
		}
	
	private File getNonCodingFileForContig(final String contig)
		{
		if(this.eigenDirectory==null) throw new IllegalStateException("Eigein directory was not define");
		int c = contig2eigen(contig);
		if(c<1) return null;
		return new File(this.eigenDirectory,"Eigen_hg19_noncoding_annot_chr"+c+".tab.bgz");
		}
	
	@Override
	public String getName()
		{
		return "Eigen";
		}

	@Override
	public String getDescription()
		{
		return "Annotator for the data of https://xioniti01.u.hpc.mssm.edu/v1.1/ : Eigen makes use of a variety of functional annotations in both coding and noncoding regions (such as made available by the ENCODE and Roadmap Epigenomics projects), and combines them into one single measure of functional importance.";
		}
	
	@Override
	public Set<VCFInfoHeaderLine> getInfoHeaderLines()
		{
		
		return new HashSet<>(Arrays.asList(
				
				));
		}

	@Override
	public Map<String, Object> getAnnotations(final VariantContext ctx)
		{
		if(!ctx.isVariant()) return Collections.emptyMap();
		final int  contig = contig2eigen(ctx.getContig());
		if( contig < 1) return Collections.emptyMap();
		final Map<String, Object> map =new HashMap<>();
		try {
			if( prev_contig==-1 || prev_contig!=contig)
				{
				CloserUtil.close(this.nonCodingFeatureReader);
				this.nonCodingFeatureReader = new TabixFeatureReader<>(
						"",new NonCodingFeatureCodec());
				prev_contig = contig;
				}
			CloseableTribbleIterator<NonCodingFeature> iter1= this.nonCodingFeatureReader.query(
						String.valueOf(contig), ctx.getStart(),ctx.getEnd());
			while(iter1.hasNext()) {
				NonCodingFeature feat = iter1.next();
				if(feat==null) continue;
				}
			iter1.close();
			
			
			CloseableTribbleIterator<CodingFeature> iter2= this.codingFeatureReader.query(
					String.valueOf(contig), ctx.getStart(),ctx.getEnd());
			while(iter2.hasNext()) {
				CodingFeature feat = iter2.next();
				if(feat==null) continue;
				}
			iter2.close();

			
			return map;
			}
		catch(IOException err) {
			throw new RuntimeException(err);
			}
		}
	@Override
	public void dispose()
		{
		prev_contig=-1;
		CloserUtil.close(this.nonCodingFeatureReader);
		CloserUtil.close(this.codingFeatureReader);
		this.nonCodingFeatureReader=null;
		this.codingFeatureReader=null;
		
		}
	}
