/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.gtf;


import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.UnaryOperator;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.TribbleException;
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
public class  GTFCodec extends AsciiFeatureCodec<GTFLine> {
	private static final String GFF_VERSION="##gff-version";
	private GTFHeaderImpl header=null;
	private UnaryOperator<String> contigNameConverter=S->S;

	public GTFCodec() {
		super(GTFLine.class);
		}
	
	/** set the contigName converter and return 'this' */
	public GTFCodec setContigNameConverter(final UnaryOperator<String> contigNameConverter) {
		this.contigNameConverter = (contigNameConverter==null?S->S:contigNameConverter);
		return this;
		}
	
	@Override
	public GTFLine decode(final String line) {
		/* non, on s'en fou a vrai dire...
		if(this.header==null) {
			throw new RuntimeIOException("header was not parsed");
		}*/
		if(StringUtil.isBlank(line) || line.startsWith("#")) return null;
		return decode(CharSplitter.TAB.split(line));
		}
	
	/** decode from array of strings */
	public GTFLine decode(String tokens[]) {
		if(tokens==null) return null;
		final String ctg = this.contigNameConverter.apply(tokens[0]);
		if(StringUtils.isBlank(ctg)) return null;
		if(!ctg.equals(tokens[0])) {
			tokens = Arrays.copyOf(tokens, tokens.length);
			tokens[0] = ctg;
			}
		return new GTFLineImpl(tokens);
		}
	
	
	public static interface GTFHeader
		{
		public boolean isGff3();
		public List<String> getLines();
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
			final String line = r.next();
			if(StringUtil.isBlank(line) || line.startsWith("#")) continue;
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

	
	private static class GTFLineImpl implements GTFLine
		{
		private final String tokens[];
		private final int start;
		private final int end;
		private Map<String,String> attributes = null;
		public GTFLineImpl(final String tokens[])
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
			return getAttributes().get(key);
			}
		
		@Override
		public Map<String, String> getAttributes() {
			if(this.attributes==null) {
				this.attributes = parseAttributes(get(8));
				}
			return this.attributes;
			}
		
		@Override
		public String toString() {
			return getLine();
			}
		}

	
    /**
     * Parse attributes field for GTF feature
     * 
     * Warning: I'm not sure it's the best algotithm to parse attributes UTF-8 escaping is not
     * handled for now Multiple values for the same key are defined with multiple K=V1; K=V2 Is it
     * the right gtf way ?
     * 
     * 
     * @param attributesString attributes field string from line in GTF file
     * @return map of keys to values for attributes of this feature
     * @throws UnsupportedEncodingException
     */
    static private Map<String,String> parseAttributes(final String attributesString) {
        final char ATTRIBUTE_DELIMITER = ';';
        final String NO_VALUE = ".";
    	if (attributesString.equals(NO_VALUE)) {
            return Collections.emptyMap();
        }
        final Map<String,String> attributes = new LinkedHashMap<>();
        final int len = attributesString.length();
        int i = 0;
        for (;;) {
            // skip whitespaces
            while (i < len && Character.isWhitespace(attributesString.charAt(i))) {
                i++;
            }
            // end of string
            if (i >= len) {
                break;
            }
            
            final StringBuilder keyBuilder = new StringBuilder();

            // consumme key
            while (i < len && !Character.isWhitespace(attributesString.charAt(i))) {
                keyBuilder.append(attributesString.charAt(i));
                i++;
            }
            // skip whitespaces
            while (i < len && Character.isWhitespace(attributesString.charAt(i))) {
                i++;
            }

            final String key = keyBuilder.toString();


            // no value
            if (i >= len) {
                //logger.warn("no value for '" + key + "' in " + attributesString);
            	attributes.put(key,"");
                break;
            }

            /* read VALUE */
            final StringBuilder valueBuilder = new StringBuilder();
            // first char of value
            char c = attributesString.charAt(i);

            if (c == ATTRIBUTE_DELIMITER) { // no value
                i++;
                attributes.put(key,"");
                //logger.warn("no value for '" + key + "' in " + attributesString);
                continue;
            }


            // quoted string
            if (c == '\"') {
                i++;
                while (i < len) {
                    c = attributesString.charAt(i);
                    ++i;
                    if (c == '\\') {
                        c = (i < len ? attributesString.charAt(i) : '\0');
                        ++i;
                        switch (c) {
                            case '"':
                                valueBuilder.append("\"");
                                break;
                            case '\'':
                                valueBuilder.append("\'");
                                break;
                            case 't':
                                valueBuilder.append("\t");
                                break;
                            case 'n':
                                valueBuilder.append("\n");
                                break;
                            default:
                                //logger.warn("unparsed value in " + attributesString);
                                break;
                        }
                    } else if (c == '\"') {
                        attributes.put(key,valueBuilder.toString());
                        break;
                    } else {
                        valueBuilder.append(c);
                    }
                } // end of while
                
                // skip whitespaces
                while (i < len && Character.isWhitespace(attributesString.charAt(i))) {
                    i++;
                }
                if (i >= len) {
                    break;
                }
                if (attributesString.charAt(i) != ATTRIBUTE_DELIMITER)
                    throw new TribbleException("expected a " + ATTRIBUTE_DELIMITER
                            + " after value " + valueBuilder);
                i++;
            } else /* not a quoted string */
            {
                while (i < len) {
                    c = attributesString.charAt(i);
                    ++i;
                    if (c == ATTRIBUTE_DELIMITER) {
                        break;
                    }
                    valueBuilder.append(c);
                }
                attributes.put(key,valueBuilder.toString());
            }
        }
        return attributes;
    	}
	
	}
