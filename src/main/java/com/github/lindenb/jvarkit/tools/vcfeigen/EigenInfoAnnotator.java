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
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.util.vcf.InfoAnnotator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class EigenInfoAnnotator implements InfoAnnotator
	{
	private static class KeyAndType
		{
		final String key;
		final VCFHeaderLineCount count;
		final VCFHeaderLineType type;
		
		KeyAndType(final String key,final VCFHeaderLineCount count,VCFHeaderLineType type)
			{
			this.key=key;
			this.count=count;
			this.type=type;
			
			}
		KeyAndType(final String key)
			{
			this(key,VCFHeaderLineCount.A,VCFHeaderLineType.Float);
			}
		}
	
	public static abstract class AbstractFeature implements Feature
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
			this.alt = Allele.create(tokens[3], false);
			this.tokens=tokens;
			}
		
		boolean accept(final VariantContext ctx) {
			int C1 = contig2eigen(ctx.getContig());
			int C2 = contig2eigen(this.contig);
			if(C1==-1 || C2==-1 || C1!=C2) return false;
			if(!ctx.getReference().equals(this.ref)) return false;
			if(!ctx.hasAllele(this.alt)) return false;
			return true;
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
		String get(int i) {
			return i>=0 && i< this.tokens.length?tokens[i]:null;
			}
		
		@Override
		public String toString() {
			return Arrays.toString(this.tokens);
			}
		}
	public  static class CodingFeature extends AbstractFeature
		{
		private CodingFeature(final String tokens[]) {
			super(tokens);
			}
		}
	public  static class NonCodingFeature extends AbstractFeature
		{
		private NonCodingFeature(final String tokens[]) {
			super(tokens);
			}
		}
	
	public static abstract class AbstractFeatureCodec<T extends AbstractFeature>  extends AsciiFeatureCodec<T>
		{
		protected final Pattern tab=Pattern.compile("[\t]");
		protected AbstractFeatureCodec(Class<T> C){
			super(C);
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
	
	/* public because http://gatkforums.broadinstitute.org/gatk/discussion/8851/ */
	public static class CodingFeatureCodec extends AbstractFeatureCodec<CodingFeature>
		{
		public CodingFeatureCodec() {
			super(CodingFeature.class);
			}
		@Override
		public CodingFeature decode(final String line)
			{
			if(line.startsWith("chr")) return null;
			return new CodingFeature(tab.split(line));
			}
		
		
		}
	/* public because http://gatkforums.broadinstitute.org/gatk/discussion/8851/ */
	public static class NonCodingFeatureCodec extends AbstractFeatureCodec<NonCodingFeature>
		{
		public NonCodingFeatureCodec() {
			super(NonCodingFeature.class);
			}
		@Override
		public NonCodingFeature decode(String line)
			{
			if(line.startsWith("chr")) return null;
			return new NonCodingFeature(tab.split(line));
			}
		
		}
	
	
	private File eigenDirectory = null;
	private final List<VCFInfoHeaderLine> noncodingheaderlines = new ArrayList<>();
	private final List<VCFInfoHeaderLine> codingheaderlines = new ArrayList<>();
	private TabixFeatureReader<NonCodingFeature, PositionalBufferedStream> nonCodingFeatureReader = null;
	private TabixFeatureReader<CodingFeature, PositionalBufferedStream> codingFeatureReader = null;
	private int prev_contig=-1;
	
	
	public EigenInfoAnnotator(final File dir) {
		this.eigenDirectory = dir;
		IOUtil.assertDirectoryIsReadable(dir);
		
		//this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"chr", VCFHeaderLineCount.A, VCFHeaderLineType.String, "chr"));
		//this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"position", VCFHeaderLineCount.A, VCFHeaderLineType.String, "position"));
		//this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"ref", VCFHeaderLineCount.A, VCFHeaderLineType.String, "ref"));
		//this.noncodingheaderlines.add(new VCFInfoHeaderLine(prefix+"alt", VCFHeaderLineCount.A, VCFHeaderLineType.String, "alt"));
		
		for(final KeyAndType key : new KeyAndType[]{
			new KeyAndType("GERP_NR"),	new KeyAndType("GERP_RS"),	new KeyAndType("PhyloPri"),	new KeyAndType("PhyloPla"),
			new KeyAndType("PhyloVer"),	new KeyAndType("PhastPri"),	new KeyAndType("PhastPla"),	new KeyAndType("PhastVer"),
			new KeyAndType("H3K4Me1"),	new KeyAndType("H3K4Me3"),	new KeyAndType("H3K27ac"),	new KeyAndType("TFBS_max"),
			new KeyAndType("TFBS_sum"),	new KeyAndType("TFBS_num"),	new KeyAndType("OCPval"),	new KeyAndType("DnaseSig"),
			new KeyAndType("DnasePval"),	new KeyAndType("FaireSig"),	new KeyAndType("FairePval"),	new KeyAndType("PolIISig"),
			new KeyAndType("PolIIPval"),	new KeyAndType("ctcfSig"),	new KeyAndType("ctcfPval"),	new KeyAndType("cmycSig"),
			new KeyAndType("cmycPval"),	new KeyAndType("Eigen-raw"),	new KeyAndType("Eigen-phred"),	new KeyAndType("Eigen-PC-raw"),
			new KeyAndType("Eigen-PC-phred")}
			)
			{
			final String prefix="EIGEN_NC_";
			final VCFInfoHeaderLine vihl = new VCFInfoHeaderLine(
					prefix+key.key.replace('-', '_'),
					key.count, 
					key.type,
					key.key+ " (Non-Coding)"
					);
			this.noncodingheaderlines.add(vihl);
			}
		
		for(final KeyAndType key : new KeyAndType[]{
				new KeyAndType("SIFT"),	new KeyAndType("PolyPhenDIV"),	new KeyAndType("PolyPhenVar"),	new KeyAndType("MA"),
				new KeyAndType("GERP_NR"),	new KeyAndType("GERP_RS"),	new KeyAndType("PhyloPri"),	new KeyAndType("PhyloPla"),
				new KeyAndType("PhyloVer"),	new KeyAndType("PhastPri"),	new KeyAndType("PhastPla"),	new KeyAndType("PhastVer"),
				new KeyAndType("Consequence",VCFHeaderLineCount.A,VCFHeaderLineType.String),
				new KeyAndType("Eigen-raw"),	new KeyAndType("Eigen-phred"),	new KeyAndType("Eigen-PC-raw"),
				new KeyAndType("Eigen-PC-phred")
				}
				)
				{
				final String prefix="EIGEN_CODING_";
				final VCFInfoHeaderLine vihl = new VCFInfoHeaderLine(
						prefix+key.key.replace('-', '_'),
						key.count, 
						key.type,
						key.key+ " (Coding)"
						);
				this.codingheaderlines.add(vihl);
				}
		
		}
	
	private static int contig2eigen(String chr) {
		if(chr.toLowerCase().startsWith("chr")) chr=chr.substring(3);
		if(chr.isEmpty() || !Character.isDigit(chr.charAt(0))) return -1;
		try {
			int c=Integer.parseInt(chr.trim());
			if(c<1 || c>22) return -1;
			return c;
			}
		catch(NumberFormatException err) {
			return -1;
			}
		}
	private String tabixPrefix="Eigen_hg19_";
	public void setTabixFilePrefix(String prefix) {
		this.tabixPrefix = prefix;
	}
	
	public String getTabixPrefix() {
		return tabixPrefix;
	}
	
	private File getNonCodingFileForContig(final int C)
		{
		if(this.eigenDirectory==null) throw new IllegalStateException("Eigein directory was not defined");
		
		return new File(this.eigenDirectory,getTabixPrefix()+"noncoding_annot_chr"+C+".tab.bgz");
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
		
		final Set<VCFInfoHeaderLine> set = new HashSet<>(this.noncodingheaderlines);
		set.addAll(this.codingheaderlines);
		return set;
		}

	@Override
	public Map<String, Object> getAnnotations(final VariantContext ctx)
		{
		if(!ctx.isVariant()) return Collections.emptyMap();
		final int  contig = contig2eigen(ctx.getContig());
		if( contig < 1) return Collections.emptyMap();
	
		try {
			if(this.codingFeatureReader == null) {
				this.codingFeatureReader = new TabixFeatureReader<>(
						new File(this.eigenDirectory,getTabixPrefix()+"coding_annot_04092016.tab.bgz").getPath(),
						new CodingFeatureCodec());
				}
			
			if( this.prev_contig==-1 || prev_contig!=contig)
				{
				CloserUtil.close(this.nonCodingFeatureReader);
				this.nonCodingFeatureReader = new TabixFeatureReader<>(
						getNonCodingFileForContig(contig).getPath(),
						new NonCodingFeatureCodec());
				this.prev_contig = contig;
				}
			final Map<Allele,NonCodingFeature> alt2nonCoding= new HashMap<>();
			CloseableTribbleIterator<NonCodingFeature> iter1= this.nonCodingFeatureReader.query(
						String.valueOf(contig), ctx.getStart(),ctx.getEnd());
			while(iter1.hasNext()) {
				NonCodingFeature feat = iter1.next();
				if(feat==null || !feat.accept(ctx)) continue;
				alt2nonCoding.put(feat.alt, feat);
				}
			iter1.close();
			
			final Map<Allele,CodingFeature> alt2coding= new HashMap<>();
			CloseableTribbleIterator<CodingFeature> iter2= this.codingFeatureReader.query(
					String.valueOf(contig), ctx.getStart(),ctx.getEnd());
			while(iter2.hasNext()) {
				CodingFeature feat = iter2.next();
				if(feat==null || !feat.accept(ctx)) continue;
				alt2coding.put(feat.alt, feat);
				}
			iter2.close();
			
			if( alt2nonCoding.isEmpty() && alt2coding.isEmpty() ) return Collections.emptyMap();
			
			final List<Allele> alternateAlleles =  ctx.getAlternateAlleles();
			final Map<String, Object> map =new HashMap<>();
			for(int side=0;side<2;++side)
				{
				for(int i=0;i< (side==0?noncodingheaderlines.size():codingheaderlines.size());++i)
					{
					final VCFInfoHeaderLine vihl = (side==0?noncodingheaderlines.get(i):codingheaderlines.get(i));
					final List<Object> atts = new ArrayList<>(alternateAlleles.size());
					boolean found_one=false;
					for(int altn=0;altn<alternateAlleles.size();++altn)
						{
						final Allele alt= alternateAlleles.get(altn);
						final AbstractFeature feat = (side==0?alt2nonCoding.get(alt):alt2coding.get(alt));
						if( feat == null) {
							atts.add(".");
							continue;
							}
						final String token = feat.get(4+i);
						if(token==null || token.isEmpty() || token.equals("."))
							{
							atts.add(".");
							continue;
							}
						else if( vihl.getType()==VCFHeaderLineType.String) {
							found_one=true;
							atts.add(token.replace(',', '|'));
							}
						else
							{
							try {
								Float dbl = new Float(token);
								atts.add(dbl);
								found_one=true;
							} catch (NumberFormatException err) {
								throw new IOException("Cannot cast "+token+" to float for "+vihl);
								}
							}
						}
					if(!found_one)
						{
						//nothing
						}
					else if(vihl.getCountType()==VCFHeaderLineCount.R)
						{
						map.put(vihl.getID(),new ArrayList<>(new LinkedHashSet<>(atts)));
						}
					else
						{
						if(atts.size()!=ctx.getAlternateAlleles().size()) {
							System.err.println("Number of eigen data!=number of ALTS");
						}
						map.put(vihl.getID(), atts);
						}
					}
				}
			
			return map;
			}
		catch(final IOException err) {
			throw new RuntimeException(err);
			}
		}
	@Override
	public void close()
		{
		prev_contig=-1;
		CloserUtil.close(this.nonCodingFeatureReader);
		CloserUtil.close(this.codingFeatureReader);
		this.nonCodingFeatureReader=null;
		this.codingFeatureReader=null;
		
		}
	}
