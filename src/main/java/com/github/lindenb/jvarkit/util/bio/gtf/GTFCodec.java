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
package com.github.lindenb.jvarkit.util.bio.gtf;


import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;

/**
 * GFF/GTF Codec
 * @author lindenb

Example: GFF3


JH561405.1	.	biological_region	974899	974901	0.999	+	.	logic_name=eponine
JH561405.1	.	biological_region	975262	976053	1	+	.	external_name=rank %3D 1;logic_name=firstef
JH561405.1	.	biological_region	985463	985464	0.999	-	.	logic_name=eponine
JH561405.1	ensembl	gene	987375	987680	.	+	.	ID=gene:ENSDNOG00000034810;biotype=protein_coding;gene_id=ENSDNOG00000034810;logic_name=ensembl;version=1
JH561405.1	ensembl	mRNA	987375	987680	.	+	.	ID=transcript:ENSDNOT00000032134;Parent=gene:ENSDNOG00000034810;biotype=protein_coding;transcript_id=ENSDNOT00000032134;version=1
JH561405.1	ensembl	exon	987375	987680	.	+	.	Parent=transcript:ENSDNOT00000032134;Name=ENSDNOE00000426775;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSDNOE00000426775;rank=1;version=1
JH561405.1	ensembl	CDS	987375	987680	.	+	0	ID=CDS:ENSDNOP00000020463;Parent=transcript:ENSDNOT00000032134;protein_id=ENSDNOP00000020463
JH561405.1	.	biological_region	990500	990572	24.4	-	.	external_name=Pseudo;logic_name=trnascan
JH561405.1	.	biological_region	994771	994841	36	-	.	external_name=Gly;logic_name=trnascan
JH561405.1	.	biological_region	996321	996392	23.8	-	.	external_name=Pseudo;logic_name=trnascan


Example: GTF 

JH568905.1	ensembl	transcript	3505651	3506703	.	+	.	gene_id "ENSDNOG00000033204"; gene_version "1"; transcript_id "ENSDNOT00000051349"; transcript_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";
JH568905.1	ensembl	exon	3505651	3506703	.	+	.	gene_id "ENSDNOG00000033204"; gene_version "1"; transcript_id "ENSDNOT00000051349"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSDNOE00000435289"; exon_version "1";
JH568905.1	ensembl	CDS	3505651	3506700	.	+	0	gene_id "ENSDNOG00000033204"; gene_version "1"; transcript_id "ENSDNOT00000051349"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSDNOP00000022248"; protein_version "1";
JH568905.1	ensembl	stop_codon	3506701	3506703	.	+	0	gene_id "ENSDNOG00000033204"; gene_version "1"; transcript_id "ENSDNOT00000051349"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";
JH568905.1	ensembl	gene	3667491	3667796	.	-	.	gene_id "ENSDNOG00000013511"; gene_version "1"; gene_source "ensembl"; gene_biotype "protein_coding";
JH568905.1	ensembl	transcript	3667491	3667796	.	-	.	gene_id "ENSDNOG00000013511"; gene_version "1"; transcript_id "ENSDNOT00000013503"; transcript_version "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";
JH568905.1	ensembl	exon	3667491	3667796	.	-	.	gene_id "ENSDNOG00000013511"; gene_version "1"; transcript_id "ENSDNOT00000013503"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSDNOE00000151550"; exon_version "1";
JH568905.1	ensembl	CDS	3667494	3667796	.	-	0	gene_id "ENSDNOG00000013511"; gene_version "1"; transcript_id "ENSDNOT00000013503"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSDNOP00000010472"; protein_version "1";
JH568905.1	ensembl	start_codon	3667794	3667796	.	-	0	gene_id "ENSDNOG00000013511"; gene_version "1"; transcript_id "ENSDNOT00000013503"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";
JH568905.1	ensembl	stop_codon	3667491	3667493	.	-	0	gene_id "ENSDNOG00000013511"; gene_version "1"; transcript_id "ENSDNOT00000013503"; transcript_version "1"; exon_number "1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding";


 */
public abstract class  GTFCodec extends AsciiFeatureCodec<GTFLine>{
	private static final Pattern tab=Pattern.compile("[\t]");
	private static final String GFF_VERSION="##gff-version";
	private GTFHeaderImpl header=null;
	private final Format format;
	public enum Format {gtf,gff3};
	
	public static class FormatChooser
		{
		
		@Parameter(names={"--input-gtf-format","--input-gff-format"},
				description="GTF/GFF format. Will (hopefully) determine how the gft/gff will be parsed",
				hidden=true
				)
		private Format format = Format.gff3;
		
		/** creates a new GTFCodec according to the specified format */
		public GTFCodec makeCodec() {
			switch(this.format) {
				case gff3: return new Gff3Codec();
				case gtf:
				default: return new GtfCodec();
				}
			}
		}
	
	public static interface GTFHeader
		{
		public boolean isGff3();
		public List<String> getLines();
		}

	protected static class GtfCodec extends GTFCodec
		{
		GtfCodec() {
			super(Format.gtf);
			}
		@Override
		public GTFLine decode(final String line) {
			/* non, on s'en fout à vrai dire...
			if(this.header==null) {
				throw new RuntimeIOException("header was not parsed");
			}*/
			if(line.startsWith("#") || line.isEmpty()) return null;
			return new GTFLineImpl(GTFCodec.tab.split(line));
			}
		
		}
	protected static class Gff3Codec extends GTFCodec
		{
		Gff3Codec() {
			super(Format.gff3);
			}
		@Override
		public GTFLine decode(final String line) {
			/* non, on s'en fout à vrai dire...
			if(this.header==null) {
				throw new RuntimeIOException("header was not parsed");
			}*/
			if(line.startsWith("#") || line.isEmpty()) return null;
			return new GFF3LineImpl(GTFCodec.tab.split(line));
			}
		}
	
	
	/** implementation of header */
	public static class GTFHeaderImpl implements GTFHeader
		{
		private final List<String> lines = new ArrayList<>();
		@Override
		public boolean isGff3() {
				for(final String line:this.lines) {
				if(line.startsWith(GFF_VERSION+" "))
					{
					final String version =line.substring(GFF_VERSION.length()).trim();
					if(version.equals("3"))
						{
						return true;
						}
					}
				}
			return false;
			}
		@Override
		public List<String> getLines() {
			return this.lines;
			}
		@Override
		public String toString() {
			return String.join("\n", this.lines);
			}
		}

	
	protected GTFCodec(final Format format) {
		super(GTFLine.class);
		this.format = format;
		}
	
	@Override
	public boolean canDecode(final String path) {
		if(StringUtil.isBlank(path)) return false;
		return true;
		}

	@Override
	public GTFLine decode(final LineIterator r)
		{
		for(;;)
			{
			if(!r.hasNext()) return null;
			final String line=r.next();
			if(line.startsWith("#")) continue;
			final GTFLine record =  decode(line);	
			if(record==null) continue;
			return record;
			}
		}
		
	@Override
	public GTFHeader readActualHeader(final LineIterator r) {
		if(this.header!=null) throw new RuntimeIOException("Reader already read");
		this.header = new GTFHeaderImpl();
		while(r.hasNext() && r.peek().startsWith("#"))
			{
			
			this.header.lines.add(r.next());
			}
		return this.header;
		}
		
		
	@Override
	public abstract GTFLine decode(final String line);
	
	
	
	private static abstract class AbstractGTFLineImpl implements GTFLine
		{
		final String tokens[];
		final int start;
		final int end;

		public AbstractGTFLineImpl(final String tokens[])
			{
			this.tokens = tokens;
			if(tokens.length<8)
				{	
				throw new JvarkitException.TokenErrors("Expected 8 columns",tokens);
				}
			this.start = Integer.parseInt(tokens[3]);
			this.end = Integer.parseInt(tokens[4]);
			}
		
		@Override
		public int hashCode() {
			return Arrays.hashCode(this.tokens);
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof GTFLine)) return false;
			return this.getLine().equals(GTFLine.class.cast(obj).getLine());
		}

		protected String get(int col) {
			return this.tokens[col];
		}
		@Override
		public String getContig() {
			return get(0);
		}
		@Override
		public String getSource() {
			return get(1);
		}
		@Override
		public String getType() {
			return get(2);
		}
		@Override
		public int getStart() {
			return start;
		}
		
		@Override
		public int getEnd() {
			return end;
		}
		@Override
		public Double getScore() {
			return get(5).equals(".")?null:Double.parseDouble(get(5));
		}
		@Override
		public char getStrand() {
			return get(6).charAt(0);
		}
	
		@Override
		public int getPhase() {
			return (get(7).equals(".")?GTFLine.NO_PHASE:Integer.parseInt(get(7)));
		}
		
		@Override
		public String getLine() {
			return String.join("\t", this.tokens);
			}
		
		
		@Override
		public String getAttribute(final String key) {
			final Iterator<Map.Entry<String, String>> iter = getAttributeIterator();
			while(iter.hasNext()) {
				final Map.Entry<String, String> kv = iter.next();
				if(kv.getKey().equals(key)) return kv.getValue();
				}
			return null;
			}
		
		@Override
		public Map<String, String> getAttributes() {
			final Map<String, String> hash = new HashMap<String, String>();
			final Iterator<Map.Entry<String, String>> iter = getAttributeIterator();
			while(iter.hasNext()) {
				final Map.Entry<String, String> kv = iter.next();
				hash.put(kv.getKey(),kv.getValue());
				}
			return hash;
			}
		
		@Override
		public String toString() {
			return getLine();
			}
		}

	private static class GTFLineImpl extends AbstractGTFLineImpl
		{
		GTFLineImpl(final String tokens[])
			{
			super(tokens);
			}
		@Override
		public Iterator<Entry<String, String>> getAttributeIterator() {
			return new GTFAttIter(get(8));
			}
		}

	private static class GFF3LineImpl extends AbstractGTFLineImpl
		{
		GFF3LineImpl(final String tokens[])
			{
			super(tokens);
			}
		@Override
		public Iterator<Entry<String, String>> getAttributeIterator() {
			return new GFF3AttIter(get(8));
			}
		}
	
	private static abstract class AbstractAttIter 
		extends AbstractIterator<Map.Entry<String, String>> 
		{
		protected final String mapStr;
		protected int k=0;
		AbstractAttIter(final String mapStr)
			{
			this.mapStr = mapStr;
			}
		protected void skipws() {
			while( this.k < this.mapStr.length() &&
				Character.isWhitespace(this.mapStr.charAt(this.k)))
				{
				++this.k;
				}
			}
		
		
		@Override
		protected Entry<String, String> advance() {
			skipws();
			if(k>=this.mapStr.length()) return null;
			for(;;)
				{
				skipws();
				if(this.k>=this.mapStr.length()) return null;
				char c= mapStr.charAt(k);
				if(c==';') { ++k; continue;}
				/* read KEY */
				final StringBuilder sbk=new StringBuilder();
				while( this.k < mapStr.length()) {
					c= mapStr.charAt(k);
					++k;
					if(c=='=' || Character.isWhitespace(c))
						{
						break;
						}
					sbk.append(c);
					}
				/* SKIP WS */
				skipws();
				/* EQUAL SIGN */
				if( this.k < mapStr.length() && mapStr.charAt(k)=='=') {
					++k;
					}
				/* SKIP WS */
				skipws();
				
				if( this.k >= mapStr.length())
					{
					if(sbk.length()==0) return null;
					return new AbstractMap.SimpleEntry<String,String>(sbk.toString(),"");
					}
				
				/* read VALUE */
				final StringBuilder sbv=new StringBuilder();
				c=(this.k < mapStr.length()?this.mapStr.charAt(this.k):'\0');
				// quoted string
				if( c == '\"')
					{
					++this.k;
					while( this.k < mapStr.length()) {
						c= mapStr.charAt(k);
						++k;
						if(c=='\\')
							{
							c=(k < mapStr.length()?mapStr.charAt(k):'\0');
							++k;
							switch(c) {
								case '"': sbv.append("\"");break;
								case '\'': sbv.append("\'");break;
								case 't': sbv.append("\t");break;
								case 'n': sbv.append("\n");break;
								default:break;
								}
							}
						else if(c=='\"')
							{
							break;
							}
						else
							{
							sbv.append(c);
							}
						}
					}
				else
					{
					while( this.k < this.mapStr.length()) {
						c= this.mapStr.charAt(k);
						++k;
						if(c==';' || Character.isWhitespace(c))
							{
							break;
							}
						sbv.append(c);
						}
					}
				final AbstractMap.SimpleEntry<String,String> entry= new AbstractMap.SimpleEntry<String,String>(sbk.toString(),sbv.toString());
				skipws();
				if( this.k < this.mapStr.length() && this.mapStr.charAt(this.k)==';')
					{
					this.k++;
					skipws();
					}
				return entry;
				}
			}
		}


	private static class GTFAttIter extends AbstractAttIter {
		GTFAttIter(final String mapStr)
			{
			super(mapStr);
			}
		}
	private static class GFF3AttIter extends AbstractAttIter {
		GFF3AttIter(final String mapStr)
			{
			super(mapStr);
			}
		}
	}
