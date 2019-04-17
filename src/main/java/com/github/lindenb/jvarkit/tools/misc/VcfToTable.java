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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.transform.stream.StreamResult;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.igv.IgvConstants;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.util.popgen.HardyWeinbergCalculation;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;

/**

BEGIN_DOC

## Example 

```
$ cat input.ped

FAM	M10475	0	0	1	1
FAM	M10478	0	0	2	0
FAM	M10500	M10475	M10478	2	1


$ curl -s "https://raw.githubusercontent.com/arq5x/gemini/master/test/test.region.vep.vcf" | java -jar dist/vcf2table.jar -H -p input.ped

 
INFO
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
| ID              | Type    | Count | Description                                                                                                                                             |
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
| AC              | Integer |       | Allele count in genotypes, for each ALT allele, in the same order as listed                                                                             |
| AF              | Float   |       | Allele Frequency, for each ALT allele, in the same order as listed                                                                                      |
| AN              | Integer | 1     | Total number of alleles in called genotypes                                                                                                             |
| BaseQRankSum    | Float   | 1     | Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities                                                                                       |
| CSQ             | String  |       | Consequence type as predicted by VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|ALLELE_NUM |
| DP              | Integer | 1     | Approximate read depth; some reads may have been filtered                                                                                               |
| DS              | Flag    | 0     | Were any of the samples downsampled?                                                                                                                    |
| Dels            | Float   | 1     | Fraction of Reads Containing Spanning Deletions                                                                                                         |
| FS              | Float   | 1     | Phred-scaled p-value using Fisher's exact test to detect strand bias                                                                                    |
| HRun            | Integer | 1     | Largest Contiguous Homopolymer Run of Variant Allele In Either Direction                                                                                |
| HaplotypeScore  | Float   | 1     | Consistency of the site with at most two segregating haplotypes                                                                                         |
| InbreedingCoeff | Float   | 1     | Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation                       |
| MQ              | Float   | 1     | RMS Mapping Quality                                                                                                                                     |
| MQ0             | Integer | 1     | Total Mapping Quality Zero Reads                                                                                                                        |
| MQRankSum       | Float   | 1     | Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities                                                                               |
| QD              | Float   | 1     | Variant Confidence/Quality by Depth                                                                                                                     |
| ReadPosRankSum  | Float   | 1     | Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias                                                                                   |
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+

FORMAT
+----+---------+-------+----------------------------------------------------------------------------------------+
| ID | Type    | Count | Description                                                                            |
+----+---------+-------+----------------------------------------------------------------------------------------+
| AD | Integer |       | Allelic depths for the ref and alt alleles in the order listed                         |
| DP | Integer | 1     | Approximate read depth (reads with MQ=255 or with bad mates are filtered)              |
| GQ | Integer | 1     | Genotype Quality                                                                       |
| GT | String  | 1     | Genotype                                                                               |
| PL | Integer |       | Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification |
+----+---------+-------+----------------------------------------------------------------------------------------+


Dict
+-----------------------+-----------+------+
| Name                  | Length    | AS   |
+-----------------------+-----------+------+
| chr1                  | 249250621 | hg19 |
(...)
| chrX                  | 155270560 | hg19 |
| chrY                  | 59373566  | hg19 |
+-----------------------+-----------+------+

Samples
+--------+---------+--------+--------+--------+------------+
| Family | Sample  | Father | Mother | Sex    | Status     |
+--------+---------+--------+--------+--------+------------+
| FAM    | M10475  |        |        | male   | affected   |
| FAM    | M10478  |        |        | female | unaffected |
| FAM    | M10500  | M10475 | M10478 | female | affected   |
| FAM    | M128215 | M10500 |        | male   | unaffected |
+--------+---------+--------+--------+--------+------------+

>>chr1/10001/T (n 1)
 Variant
 +--------+--------------------+
 | Key    | Value              |
 +--------+--------------------+
 | CHROM  | chr1               |
 | POS    | 10001              |
 | end    | 10001              |
 | ID     | .                  |
 | REF    | T                  |
 | ALT    | TC                 |
 | QUAL   | 175.91000000000003 |
 | FILTER |                    |
 | Type   | INDEL              |
 +--------+--------------------+
 Alleles
 +-----+-----+-----+-------+--------+----+----+-----+-------------+---------------+---------+-----------+
 | Idx | REF | Sym | Bases | Length | AC | AN | AF  | AC_affected | AC_unaffected | AC_male | AC_female |
 +-----+-----+-----+-------+--------+----+----+-----+-------------+---------------+---------+-----------+
 | 0   | *   |     | T     | 1      | 4  | 8  | 0.5 | 2           | 1             | 1       | 2         |
 | 1   |     |     | TC    | 2      | 4  | 8  | 0.5 | 2           | 1             | 1       | 2         |
 +-----+-----+-----+-------+--------+----+----+-----+-------------+---------------+---------+-----------+
 INFO
 +----------------+-------+----------+
 | key            | Index | Value    |
 +----------------+-------+----------+
 | AC             |       | 4        |
 | AF             |       | 0.50     |
 | AN             |       | 8        |
 | BaseQRankSum   |       | 4.975    |
 | DP             |       | 76       |
 | FS             |       | 12.516   |
 | HRun           |       | 0        |
 | HaplotypeScore |       | 218.6157 |
 | MQ             |       | 35.31    |
 | MQ0            |       | 0        |
 | MQRankSum      |       | -0.238   |
 | QD             |       | 2.31     |
 | ReadPosRankSum |       | 2.910    |
 +----------------+-------+----------+
 VEP
 +--------------------------+------+----------------+------------+-----------------+--------+------------------+-----------------------------------------------+-------------+---------+-----------------+----------------------+
 | PolyPhen                 | EXON | SIFT           | ALLELE_NUM | Gene            | SYMBOL | Protein_position | Consequence                                   | Amino_acids | Codons  | Feature         | BIOTYPE              |
 +--------------------------+------+----------------+------------+-----------------+--------+------------------+-----------------------------------------------+-------------+---------+-----------------+----------------------+
 | probably_damaging(0.956) | 8/9  | deleterious(0) | 1          | ENSG00000102967 | DHODH  | 346/395          | missense_variant                              | R/W         | Cgg/Tgg | ENST00000219240 | protein_coding       |
 |                          | 3/4  |                | 1          | ENSG00000102967 | DHODH  |                  | non_coding_exon_variant&nc_transcript_variant |             |         | ENST00000571392 | retained_intron      |
 |                          |      |                | 1          | ENSG00000102967 | DHODH  |                  | downstream_gene_variant                       |             |         | ENST00000572003 | retained_intron      |
 |                          |      |                | 1          | ENSG00000102967 | DHODH  |                  | downstream_gene_variant                       |             |         | ENST00000573843 | retained_intron      |
 |                          |      |                | 1          | ENSG00000102967 | DHODH  |                  | downstream_gene_variant                       |             |         | ENST00000573922 | processed_transcript |
 |                          |      |                | 1          | ENSG00000102967 | DHODH  | -/193            | intron_variant                                |             |         | ENST00000574309 | protein_coding       |
 | probably_damaging(0.946) | 8/9  | deleterious(0) | 1          | ENSG00000102967 | DHODH  | 344/393          | missense_variant                              | R/W         | Cgg/Tgg | ENST00000572887 | protein_coding       |
 +--------------------------+------+----------------+------------+-----------------+--------+------------------+-----------------------------------------------+-------------+---------+-----------------+----------------------+
 Genotypes
 +---------+------+-------+----+----+-----+---------+
 | Sample  | Type | AD    | DP | GQ | GT  | PL      |
 +---------+------+-------+----+----+-----+---------+
 | M10475  | HET  | 10,2  | 15 | 10 | 0/1 | 25,0,10 |
 | M10478  | HET  | 10,4  | 16 | 5  | 0/1 | 40,0,5  |
 | M10500  | HET  | 10,10 | 21 | 7  | 0/1 | 111,0,7 |
 | M128215 | HET  | 15,5  | 24 | 0  | 0/1 | 49,0,0  |
 +---------+------+-------+----+----+-----+---------+
 TRIOS
 +-----------+-----------+-----------+-----------+----------+----------+-----------+
 | Father-ID | Father-GT | Mother-ID | Mother-GT | Child-ID | Child-GT | Incompat. |
 +-----------+-----------+-----------+-----------+----------+----------+-----------+
 | M10475    | 0/1       | M10478    | 0/1       | M10500   | 0/1      |           |
 +-----------+-----------+-----------+-----------+----------+----------+-----------+
<<chr1/10001/T n 1

(...)
```

### Html output:


```
$ java -jar dist/vcf2table.jar file.vcf --color --format html > out.html
```

### Screenshots

[https://twitter.com/yokofakun/status/1067730485487366145](https://twitter.com/yokofakun/status/1067730485487366145)

![https://pbs.twimg.com/media/DtFXEhLWkAE5roc.jpg](https://pbs.twimg.com/media/DtFXEhLWkAE5roc.jpg)

[https://twitter.com/yokofakun/status/922475502933368832](https://twitter.com/yokofakun/status/922475502933368832)

![https://pbs.twimg.com/media/DM1KdWFX0AUfbxR.jpg](https://pbs.twimg.com/media/DM1KdWFX0AUfbxR.jpg)

END_DOC

## History

* 20190416 : colors for AD
* 20170517 : fixed a bug in removeEmptyColumns

 */
@Program(name="vcf2table",
		description="convert a vcf to a table, to ease display in the terminal",
		keywords={"vcf","table","visualization"},
		references="Vcf2table : a VCF prettifier. Lindenbaum & al. 2018. figshare. [https://doi.org/10.6084/m9.figshare.5853801](https://doi.org/10.6084/m9.figshare.5853801)",
		biostars=293855,
		modificationDate="20190416"
		)
public class VcfToTable extends Launcher {
	private static final Logger LOG = Logger.build(VcfToTable.class).make();
	private static final String DEFAULT_MARGIN=" ";
	public static final String ANSI_ESCAPE = "\u001B[";
	public static final String ANSI_RESET = ANSI_ESCAPE+"0m";
	public static final String DEFAULT_CSS_STYLE="div.variant0 {background-color:#E4FAFF;}\ndiv.variant1 {background-color:#C6DDE2;}\ntable.minimalistBlack { border: 1px solid #1F1F1F; text-align: left; border-collapse: collapse; } table.minimalistBlack td, table.minimalistBlack th { border: 1px solid #1F1F1F; padding: 5px 2px; } table.minimalistBlack tbody td { font-size: 13px; } table.minimalistBlack thead { background: #CFCFCF; background: -moz-linear-gradient(top, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%); background: -webkit-linear-gradient(top, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%); background: linear-gradient(to bottom, #dbdbdb 0%, #d3d3d3 66%, #CFCFCF 100%); border-bottom: 2px solid #000000; } table.minimalistBlack thead th { font-size: 15px; font-weight: bold; color: #000000; text-align: left; } table.minimalistBlack tfoot td { font-size: 14px; } ";

	/** public: can be used by tools using vcf2table */
	public enum OutputFormat {
		text,html
		}
	
	private enum AnsiColor {
    	BLACK (30),
    	RED (31),
    	GREEN (32),
    	YELLOW (33),
    	BLUE (34),
    	MAGENTA (35),
    	CYAN (36),
    	WHITE (37)
		;
    	
    	AnsiColor(final int opcode) {
    		this.opcode=opcode;
    		}
    	final int opcode;
    	}

	
	private static abstract class Decorator
		{
		final Object o;
		Decorator(final Object o)
			{
			this.o = o;
			}
		
		protected String escapeHttp(final String s) {
			try {
				return StringUtils.escapeHttp(s);
				}
			catch(final Exception err)
				{
				return s;
				}
			}
		
		void beginXml(XMLStreamWriter w) throws XMLStreamException  {}
		void endXml(XMLStreamWriter w) throws XMLStreamException  {}
		
		void writeXml(XMLStreamWriter w) throws XMLStreamException
			{
			if(o==null) {
				return;
				}
			else if(o instanceof Decorator)
				{
				final Decorator xo=Decorator.class.cast(o);
				beginXml(w);
				xo.writeXml(w);
				endXml(w);
				}
			else
				{
				final String s= this.toString();
				if(!StringUtil.isBlank(s)) {
					beginXml(w);
					w.writeCharacters(s);
					endXml(w);
					}
				}
			}
		public String start() { return "";}
		public String end() { return "";}
		void print(final PrintStream out) {
			if(o==null) {
				return;
				}
			else if(o instanceof Decorator)
				{
				final Decorator xo=Decorator.class.cast(o);
				out.print(start());
				xo.print(out);
				out.print(end());
				}
			else
				{
				final String s= this.toString();
				if(!StringUtil.isBlank(s)) {
					out.print(start());
					out.print(o.toString());
					out.print(end());
					}
				}
			}
		
		@Override
		public String toString() {
			return o==null?"":o.toString();
			}
		}
	
	private static class HyperlinkDecorator extends Decorator
		{
		HyperlinkDecorator(final Object o) {
			super(o);
			}
		
		protected String getURL()  	{
			final String str= this.toString();
			if(StringUtil.isBlank(str)) return null;
			if(Pattern.compile("rs[0-9]+").matcher(str.toLowerCase()).matches())
				{
				return "https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs="+str.substring(2);
				}
			else if(
					Pattern.compile("ENS[TPGR][0-9]+").matcher(str.toUpperCase()).matches() ||
					Pattern.compile("ENSEST[TGP][0-9]+").matcher(str.toUpperCase()).matches() 

					)
				{
				return "http://www.ensembl.org/Multi/Search/Results?species=all;idx=;q="+str.toUpperCase()+";species=;site=ensembl";
				}
			else if(Pattern.compile("[XN][MR]_[0-9\\.]+").matcher(str.toUpperCase()).matches())
				{
				return "https://www.ncbi.nlm.nih.gov/nuccore/"+str;
				}
			else if(Pattern.compile("[NX]P_[0-9\\.]+").matcher(str.toUpperCase()).matches())
				{
				return "https://www.ncbi.nlm.nih.gov/protein/"+str;
				}
			else if(Pattern.compile("CCDS[0-9\\.]+").matcher(str.toUpperCase()).matches())
				{
				return "https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&GO=MainBrowse&DATA="+str;
				}
			else if(IOUtil.isUrl(str))
				{
				return str;
				}
			return null;
			}
		
		
		
		@Override
		void beginXml(final XMLStreamWriter w) throws XMLStreamException {
			final String u = this.getURL();
			if(StringUtil.isBlank(u)) return;
			w.writeStartElement("a");
			w.writeAttribute("title",this.toString());
			w.writeAttribute("target","_blank");
			w.writeAttribute("href",u);
			}
		@Override
		void endXml(final XMLStreamWriter w) throws XMLStreamException {
			final String u = this.getURL();
			if(StringUtil.isBlank(u)) return;
			w.writeEndElement();
			}
		}
	
	
	private static class SODecorator extends Decorator
		{
		SODecorator(final String s)
			{
			super(s);
			}
		@Override
		void writeXml(XMLStreamWriter w) throws XMLStreamException {
			final String str= this.toString();
			if(StringUtil.isBlank(str)) return;
			final String tokens[]=str.split("[;&]");
			for(int i=0;i< tokens.length;i++)
				{
				if(i>0) w.writeCharacters("; ");
				final String s=tokens[i];
				if(StringUtil.isBlank(s)) continue;
				w.writeStartElement("a");
				w.writeAttribute("href","http://www.sequenceontology.org/browser/obob.cgi?release=current_svn&rm=term_list&obo_query="+escapeHttp(s));
				w.writeAttribute("target","_blank");
				w.writeAttribute("title",s);
				w.writeCharacters(s);
				w.writeEndElement();//a
				}
			}
		}
	
	private static class GenelinkDecorator extends HyperlinkDecorator
		{
		GenelinkDecorator(final String gene) {
			super(gene);
			}
		protected String getURL() {
			String u= super.getURL();
			if(!StringUtil.isBlank(u)) return u;
			final String str= this.toString();
			if(StringUtil.isBlank(str)) return null;
			return "https://www.ncbi.nlm.nih.gov/gene/?term="+escapeHttp(str);
			}
		}
	
	private static class HgncDecorator extends HyperlinkDecorator
		{
		HgncDecorator(final String gene) {
			super(gene);
			}
		protected String getURL() {
			final String str= this.toString();
			if(StringUtil.isBlank(str)) return null;
			return "https://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id="+escapeHttp(str);
			}
		}

	
	private static class ColoredDecorator extends Decorator
		{
		final AnsiColor color;
		ColoredDecorator(final Object o,final AnsiColor color) {super(o);this.color=color;}
		@Override public String start() {
			if(o==null || color==null) return "";
			return ANSI_ESCAPE+this.color.opcode+"m";
			}
		@Override public String end()
			{
			if(o==null || color==null) return "";
			return ANSI_RESET;
			}
		@Override
		void beginXml(final XMLStreamWriter w) throws XMLStreamException {
			if(o==null || color==null) return;
			w.writeStartElement("span");
			w.writeAttribute("style", "color:"+this.color.name().toLowerCase());
			}
		@Override
		void endXml(final XMLStreamWriter w) throws XMLStreamException {
			if(o==null || color==null) return;
			w.writeEndElement();
			}
		}

	
	private static class Column
		{
		final String label;
		int maxLength=0;
		Column(final String label) {
			this.label=label;
			this.maxLength=label.length();
			}
		}

	
	
	private static class Table
		{
		private String caption="";
		private final List<Column> columns;
		private final List<List<Object>> rows = new ArrayList<>();
		public Table(final String...header) {
			this(Arrays.asList(header));
			}
		public Table(final List<String> labels)
			{
			this.columns= new ArrayList<>(labels.stream().
					map(L->new Column(L)).
					collect(Collectors.toList()));
			}
		public Table setCaption(final String t) {
			this.caption=t;
			return this;
		}
		public void addRow(final Object...items)
			{	
			addList(Arrays.asList(items));
			}
		public void addList(final List<Object> row)
			{
			this.rows.add(new ArrayList<>(row));
			for(int i=0;i< row.size()&& i< this.columns.size() ;++i)
				{
				final Object o = row.get(i);
				if(o==null) continue;
				final String str=o.toString();
				final int len = str.length();
				this.columns.get(i).maxLength = Math.max(this.columns.get(i).maxLength,len);
				}
			}
		
		public Table removeEmptyColumns()
			{
			int i=0;
			while(i < this.columns.size())
				{
				
				boolean empty=true;
				for(final List<Object> row:this.rows)
					{
					if(i>=row.size()) continue;
					final Object o= row.get(i);
					if(o==null || o.equals("")) continue;
					empty=false;
					break;
					}
				if(empty)
					{
					this.columns.remove(i);
					for(final List<Object> row:this.rows)
						{
						if(i>=row.size()) continue;
						row.remove(i);
						}
					}
				else
					{
					++i;
					}
				}
			return this;
			}
		/*
		private int tableWidth() {
			return this.columns.stream().mapToInt(C->C.maxLength +2 ).sum() 
					+ 1
					+ this.columns.size();
			}*/
	
		
		public void print(final String margin,final PrintStream out) {
			if(this.rows.isEmpty() || this.columns.isEmpty()) return;
			final StringBuilder hr= new StringBuilder();
			
			out.print(margin);
			
			out.println(this.caption);
			hr.append('+');
			for(int i=0;i< this.columns.size();i++)
				{
				hr.append("-");
				for(int j=0;j< this.columns.get(i).maxLength;j++)hr.append('-');
				hr.append("-+");
				}
			out.print(margin);
			out.println(hr.toString());
			
			out.print(margin);
			out.print('|');
			for(int i=0;i< this.columns.size();i++)
				{
				out.print(" ");
				out.print(this.columns.get(i).label);
				for(int j=this.columns.get(i).label.length();j< this.columns.get(i).maxLength;j++) out.print(' ');
				out.print(" |");
				}
			out.println();
			
			out.print(margin);
			out.println(hr.toString());
			
			for(int y=0;y< this.rows.size();++y) {
				final List<Object> row= this.rows.get(y);
				out.print(margin);
				out.print("|");
				for(int i=0;i< this.columns.size() && i< row.size();i++)
					{
					final Object cell = row.get(i);
					final String str=  cell==null?"":cell.toString();
					out.print(" ");
					if( cell instanceof Decorator)
						{
						Decorator.class.cast(cell).print(out);
						}
					else
						{
						out.print(str);
						}
					for(int j=str.length();j< this.columns.get(i).maxLength;j++) out.print(' ');
					out.print(" |");
					}
				out.println();
				
			}
			out.print(margin);
			out.println(hr.toString());
			}
		
		void writeText(final XMLStreamWriter w,final Object o) throws XMLStreamException
			{
			if(o==null) return ;
			w.writeCharacters(String.valueOf(o));
			}
		
		void write(final XMLStreamWriter w) throws XMLStreamException
			{
			if(this.rows.isEmpty() || this.columns.isEmpty()) return;
			w.writeStartElement("div");
			

			w.writeStartElement("table");
			w.writeAttribute("class","minimalistBlack");
			w.writeStartElement("thead");
			w.writeStartElement("caption");
			writeText(w,this.caption);
			w.writeEndElement();//caption
			w.writeStartElement("tr");
			for(int i=0;i< this.columns.size();i++)
				{
				w.writeStartElement("td");
				writeText(w,this.columns.get(i).label);
				w.writeEndElement();//td
				}
			w.writeEndElement();//tr
			w.writeEndElement();//thead
			w.writeStartElement("tbody");
			
			
			
			
			for(int y=0;y< this.rows.size();++y) {
				w.writeStartElement("tr");
				final List<Object> row= this.rows.get(y);
				
				for(int i=0;i< this.columns.size() && i< row.size();i++)
					{
					w.writeStartElement("td");
					final Object cell = row.get(i);
					final String str=  cell==null?"":cell.toString();
					
					if(cell instanceof Decorator) 
						{
						Decorator.class.cast(cell).writeXml(w);
						}
					else
						{
						writeText(w,str);
						}
					w.writeEndElement();//td
					}
				w.writeEndElement();//tr
				
				}
			w.writeEndElement();//tbody
			w.writeEndElement();//table
			w.writeEndElement();//div
			}
		
		}
	
	
	
	public static class VcfToTableViewer
		implements VariantContextWriter
		{
		@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
		private File outputFile = null;
		@Parameter(names={"-H","--header"},description="Print Header")
		private boolean printHeader=false;
		@Parameter(names={"-g","--hideGenotypes"},description="Hide All genotypes")
		private boolean hideGenotypes=false;
		@Parameter(names={"-nc","--hideNoCalls"},description="Hide NO_CALL genotypes")
		private boolean hideNoCallGenotypes=false;
		@Parameter(names={"-hr","--hideHomRefs"},description="Hide HOM_REF genotypes")
		private boolean hideHomRefGenotypes=false;
		@Parameter(names={"-p","--ped","--pedigree"},description="Optional Pedigree file:"+Pedigree.OPT_DESCRIPTION+" If undefined, this tool will try to get the pedigree from the header.")
		private File pedigreeFile=null;
		@Parameter(names={"-L","-limit","--limit"},description="Limit the number of output variant. '-1' == ALL/No limit.")
		private int limitVariants=-1;
		@Parameter(names={"--color","--colors"},description="[20170808] Print Terminal ANSI colors.")
		private boolean useANSIColors=false;
		@Parameter(names={"--hideAlleles"},description="[20170808] hide Alleles table.")
		private boolean hideAlleles=false;
		@Parameter(names={"--hideFilters"},description="[20170808] hide Filters table.")
		private boolean hideFilters=false;
		@Parameter(names={"--hideInfo"},description="[20170808] hide INFO table.")
		private boolean hideInfo=false;
		@Parameter(names={"--hidePredictions"},description="[20170808] hide SNPEFF/VEP table.")
		private boolean hidePredictions=false;
		@Parameter(names={"--hideGTypes"},description="[20180221] hide Genotype.Type table")
		private boolean hideGTypes=false;
		@Parameter(names={"--hideHyperlinks"},description="[20191102] hide Hyperlinks table.")
		private boolean hideHyperlinks=false;

		@Parameter(names={"--format"},description="[20171020] output format.")
		private OutputFormat outputFormat=OutputFormat.text;
		@Parameter(names={"--no-html-header"},description="[20171023] ignore html header for HTML output.")
		private boolean hideHtmlHeader=false;
		@Parameter(names={"--url"},description=Launcher.USER_CUSTOM_INTERVAL_URL_DESC)
		private String userCustomUrl=null;
		@Parameter(names={"--google"},description="use google charts (HTML only)")
		private boolean googleChart  = false;
		@Parameter(names={"--chartsize"},description="google charts dimension (HTML only). Format (integer)x(integer). eg: '1000x500' or (width) e.g: '1000'")
		private String googleChartSizeStr  = null;

		
		private AbstractViewer delegate=null;
		private PrintStream outputStream = System.out;
		
		private abstract class AbstractViewer implements VariantContextWriter
			{
			protected int countVariants=0;
			private VCFHeader header=null;
			private VCFEncoder vcfEncoder=null;
			private Pedigree pedigree = null;
			private VcfTools vcfTools = null;

			VcfToTableViewer getOwner() {
				return VcfToTableViewer.this;
				}
			
			abstract void writeTable(final String margin,final Table t);
		
			abstract void println(String s);
			abstract void println();
			abstract void startVariant(final VariantContext ctx);
			abstract void endVariant(final VariantContext ctx);
			
			
			@Override
			public void writeHeader(final VCFHeader header)
				{
				this.header = header;
				this.vcfEncoder = new VCFEncoder(header, true, true);
				this.vcfTools = new VcfTools(header);
				
				
				if(getOwner().pedigreeFile!=null) {
					try {
						this.pedigree = Pedigree.newParser().parse(getOwner().pedigreeFile);
					} catch (final IOException e) {
						throw new RuntimeIOException(e);
						}
				} else
					{
					this.pedigree = Pedigree.newParser().parse(header);
					}
				
				if(this.pedigree.isEmpty())
					{
					this.pedigree = null;
					}
				String margin="";
				
				
				if(getOwner().printHeader)
					{
					/** INFO */
						{
						final Table t=new Table("ID","Type","Count","CountType","Description").setCaption("INFO");
						header.getInfoHeaderLines().stream().
							map(F->{
									final List<Object> r=new ArrayList<>(5);
									r.add(F.getID());
									r.add(F.getType()==null?null:F.getType().name());
									r.add(F.isFixedCount()?F.getCount():null);
									r.add(F.getCountType());
									r.add(F.getDescription());
									return r;
									}).
							forEach(R->t.addList(R));
						writeTable(margin,t);
						this.println();
						}
					/** FORMAT */
						{
						final Table t=new Table("ID","Type","Count","CountType","Description").setCaption("FORMAT");
						header.getFormatHeaderLines().stream().
							map(F->{
									final List<Object> r=new ArrayList<>(5);
									r.add(F.getID());
									r.add(F.getType()==null?null:F.getType().name());
									r.add(F.isFixedCount()?F.getCount():null);
									r.add(F.getCountType());
									r.add(F.getDescription());
									return r;
									}).
							forEach(R->t.addList(R));
						this.writeTable(margin,t);
						this.println();
						}
					/** FILTER */
						{
						final Table t=new Table("ID","Description").setCaption("FILTERS");
						header.getFilterLines().forEach(
							L->t.addRow(L.getID(),L.getDescription())
							);
						this.writeTable(margin,t.removeEmptyColumns());
						this.println();
						}
					
					/** OTHER METADATA */
						{
						final Table t=new Table("ID","Description").setCaption("Metadata");
						header.getOtherHeaderLines().forEach(
							L->t.addRow(L.getKey(),L.getValue())
							);
						
						this.writeTable(margin,t.removeEmptyColumns());
						this.println();
						}
						
					/** DICT */
						{
						final SAMSequenceDictionary dict = header.getSequenceDictionary();
						if (dict != null) {
							final List<String> h = new ArrayList<>();
							h.add("Name");
							h.add("Length");
							final Set<String> all_attributes = dict.getSequences().stream()
									.flatMap(S -> S.getAttributes().stream()).map(A -> A.getKey())
									.filter(S -> !(S.equals("Name") || S.equals("length") || S.equals("Length"))).collect(Collectors.toSet());
							h.addAll(all_attributes);
							final Table t2 = new Table(h).setCaption("Dict");
		
							for (final SAMSequenceRecord ssr : dict.getSequences()) {
								final List<Object> r = new ArrayList<>();
								r.add(ssr.getSequenceName());
								r.add(ssr.getSequenceLength());
								for (final String key : all_attributes) {
									r.add(ssr.getAttribute(key));
								}
								t2.addList(r);
							}
							this.writeTable(margin,t2);
							}
						}
					if(this.pedigree!=null)
						{
						final Table t=new Table("Family","Sample","Father","Mother","Sex","Status").setCaption("Samples");
						for(final String sample: this.header.getSampleNamesInOrder())
							{
							final List<Object> r = new ArrayList<>();
							final Pedigree.Person person = this.pedigree.getPersonById(sample);
							
							r.add(person==null?null:person.getFamily().getId());
							r.add(sample);
							r.add(person==null?null:person.getFatherId());
							r.add(person==null?null:person.getMotherId());
							r.add(person==null?null:person.getSex());
							r.add(person==null?null:person.getStatus());
							t.addList(r);
							}
						this.writeTable(margin,t.removeEmptyColumns());
						}
						
						
					this.println();
					}
				}
			
			@Override
			public void add(final VariantContext vc)
				{

				if(getOwner().limitVariants!=-1)
					{
					if(this.countVariants== getOwner().limitVariants)
						{
						this.println();
						this.println("... More variants exist but they've been omitted...( limit reached)");
						return;
						}
					else if(this.countVariants> getOwner().limitVariants)
						{
						return;
						}
					}
				
				
				++countVariants;
				
				startVariant(vc);
				
				String margin=DEFAULT_MARGIN;
				{
				final Table t=new Table("Key","Value").setCaption("Variant");
				t.addRow("CHROM",vc.getContig());
				t.addRow("POS",vc.getStart());
				t.addRow("end",vc.getEnd());
				t.addRow("ID",vc.hasID()?new HyperlinkDecorator(vc.getID()):".");
				t.addRow("REF",vc.getReference().getDisplayString());
				t.addRow("ALT",vc.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(",")));
				t.addRow("QUAL",vc.hasLog10PError()?vc.getPhredScaledQual():null);
				if(hideFilters) /* print filters here if 'hide filters table has been set */ {
					final String filterStr=vc.isFiltered()?vc.getFilters().stream().collect(Collectors.joining(";")):null;
					t.addRow("FILTER",filterStr!=null && useANSIColors?new ColoredDecorator(filterStr, AnsiColor.RED):filterStr);
					}
				t.addRow("Type",vc.getType());
				
				
				this.writeTable(margin, t);
				}
			
			if(!getOwner().hideAlleles)
				{
				boolean has_affected_cols=false;
				int AN=-1;
				final List<String> h = new ArrayList<>(Arrays.asList("Idx","REF","Sym","Bases","Length"));
				if(vc.hasGenotypes())
					{
					h.add("HW");
					h.add("AC");
					h.add("AN");
					h.add("AF");
					AN = (int)vc.getGenotypes().stream().flatMap(G->G.getAlleles().stream()).filter(A->A.isCalled()).count();
					
					if(this.pedigree!=null &&
						this.pedigree.getPersons().stream().filter(P->P.getStatus()!=Pedigree.Status.missing).findAny().isPresent()
						) {
						has_affected_cols=true;
						h.add("AC_affected");
						h.add("AC_unaffected");
						}
					if(this.pedigree!=null ) {
						h.add("AC_male");
						h.add("AC_female");
						}
					
					}
				
				
				 final Table t=new Table(h).
						 setCaption("Alleles");
				 for(final Allele a: vc.getAlleles())
				 	{
					final ArrayList<Object> r = new ArrayList<>(Arrays.asList(vc.getAlleleIndex(a),
							a.isReference()?"*":"",
							a.isSymbolic()?"*":"",
							a.getDisplayString(),
							a.isSymbolic()?null:a.length()
							));
					if(vc.hasGenotypes())
						{
						Double hw =null;
						if(!(a.isReference() || a.isNoCall()))
							{
							final Genotype aa = new GenotypeBuilder("dummy", Arrays.asList(vc.getReference(),vc.getReference())).make();
							final Genotype ab = new GenotypeBuilder("dummy", Arrays.asList(vc.getReference(),a)).make();
							final Genotype bb = new GenotypeBuilder("dummy", Arrays.asList(a,a)).make();
							final int obsaa= (int)(vc.getGenotypes().stream().filter(G->G.sameGenotype(aa, true)).count());
							final int obsab= (int)(vc.getGenotypes().stream().filter(G->G.sameGenotype(ab, true)).count());
							final int obsbb= (int)(vc.getGenotypes().stream().filter(G->G.sameGenotype(bb, true)).count());
							if( obsaa + obsab + obsbb >0) 
								{
								hw=HardyWeinbergCalculation.hwCalculate(obsaa,obsab,obsbb);
								if(hw<0) hw=null;
								}
						
							}
						
						r.add(hw);
						
						final int AC = (int)vc.getGenotypes().stream().flatMap(G->G.getAlleles().stream()).filter(A->A.equals(a, false)).count();
						r.add(AC);
						r.add(AN);
						r.add(AN<=0?".":String.valueOf(AC/(double)AN));
						if(has_affected_cols)
							{
							int AC_aff=  (int)vc.getGenotypes().stream().filter(G->{
										final Pedigree.Person p=this.pedigree.getPersonById(G.getSampleName());
										if(p==null || !p.isAffected()) return false;
										return true;
										}).
									flatMap(G->G.getAlleles().stream()).filter(A->A.equals(a, false)).count();
							int AC_unaff=  (int)vc.getGenotypes().stream().filter(G->{
								final Pedigree.Person p=this.pedigree.getPersonById(G.getSampleName());
								if(p==null || !p.isUnaffected()) return false;
								return true;
								}).flatMap(G->G.getAlleles().stream()).filter(A->A.equals(a, false)).count();
							r.add(AC_aff);
							r.add(AC_unaff);
							}
						if(this.pedigree!=null ) {
							int AC_male=  (int)vc.getGenotypes().stream().filter(G->{
								final Pedigree.Person p=this.pedigree.getPersonById(G.getSampleName());
								if(p==null || !p.isMale()) return false;
								return true;
								}).flatMap(G->G.getAlleles().stream()).filter(A->A.equals(a, false)).count();
							int AC_female=  (int)vc.getGenotypes().stream().filter(G->{
								final Pedigree.Person p=this.pedigree.getPersonById(G.getSampleName());
								if(p==null || !p.isFemale()) return false;
								return true;
								}).flatMap(G->G.getAlleles().stream()).filter(A->A.equals(a, false)).count();
							r.add(AC_male);
							r.add(AC_female);
							}
						}
					t.addList(r);
				 	}
				this.writeTable(margin, t);
				}
			
			if(!getOwner().hideFilters)
				{
				/* FILTER */
				final	Table t=new Table("Filter").setCaption("FILTERS");
				 for(final String f:vc.getFilters())
				 	{
					t.addRow(getOwner().useANSIColors?new ColoredDecorator(f, AnsiColor.RED):f);
				 	}
				 this.writeTable(margin, t);
				}
			
			if(!getOwner().hideHyperlinks)
				{
				final Function<VariantContext, String> ucscContig = V->V.getContig().startsWith("chr")?V.getContig():"chr"+V.getContig();
				final Function<VariantContext, String> ensemblContig = V->
						V.getContig().startsWith("chr")?
						V.getContig().substring(3):
						V.getContig()
						;
				
				
				final	Table t=new Table("Name","URL").setCaption("Hyperlinks");
				if(vc.hasID() && vc.getID().matches("rs[0-9]+"))
					{
					t.addRow("dbSNP",new HyperlinkDecorator("https://www.ncbi.nlm.nih.gov/snp/"+vc.getID()));
					t.addRow("OpenSNP",new HyperlinkDecorator("https://opensnp.org/snps/"+vc.getID()));
					if(SequenceDictionaryUtils.isHuman(header)) {
						t.addRow("clinvar",
								new HyperlinkDecorator("https://www.ncbi.nlm.nih.gov/clinvar?term="+vc.getID()+"%5BVariant%20ID%5D"));
						}
					}
				t.addRow("IGV",new HyperlinkDecorator("https://"+ IgvConstants.DEFAULT_HOST +":"+IgvConstants.DEFAULT_PORT+"/goto?locus="+
						vc.getContig()+"%3A"+vc.getStart() +"-"+vc.getEnd()
						));
			
				
				for(final String build: new String[] {"hg19","hg38"}) {
					if(build.equals("hg19") && !SequenceDictionaryUtils.isGRCh37(header)) continue;
					if(build.equals("hg38") && !SequenceDictionaryUtils.isGRCh38(header)) continue;
					
					t.addRow("UCSC "+build,new HyperlinkDecorator("http://genome.ucsc.edu/cgi-bin/hgTracks?db="+build+"&highlight="+build+"."+
						ucscContig.apply(vc) +
						"%3A"+vc.getStart() +"-"+vc.getEnd() + "&position=" +
						ucscContig.apply(vc) +
						"%3A"+ Math.max(1,vc.getStart()-50) +"-"+(vc.getEnd()+50)
						));
					}
				
				//beacon , //varsome
				for(int side=0;side<2;++side)
					{
					if(side==0 && !SequenceDictionaryUtils.isGRCh37(header)) continue;
					if(side==1 && !SequenceDictionaryUtils.isGRCh38(header)) continue;
					for(final Allele alt: vc.getAlternateAlleles())
						{
						if(vc.getReference().isSymbolic() || alt.isSymbolic()) continue;
						//https://beacon-network.org/#/search?pos=114267128&chrom=4&allele=A&ref=G&rs=GRCh37
						t.addRow("Beacon",new HyperlinkDecorator("https://beacon-network.org/#/search?chrom="+
								ensemblContig.apply(vc) +
								"&pos="+vc.getStart()+
								"&ref="+ vc.getReference().getDisplayString()+
								"&allele="+ alt.getDisplayString()+
								"&rs="+ (side==0?"GRCh37":"GRCh38")
								));
						t.addRow("Varsome",new HyperlinkDecorator("https://varsome.com/variant/"+
								(side==0?"hg19/":"hg38/")+
								ensemblContig.apply(vc) + "-"+
								vc.getStart()+ "-"+
								vc.getReference().getDisplayString()+"-"+
								alt.getDisplayString()
								));
						}
					}
				
								
				if(SequenceDictionaryUtils.isGRCh37(header)) {
					
					for(final Allele alt: vc.getAlternateAlleles())
						{
						if(vc.getReference().isSymbolic() || alt.isSymbolic()) continue;
						// marvel https://twitter.com/julawang/status/1094666160711323649
						t.addRow("Marrvel",new HyperlinkDecorator("http://marrvel.org/search/variant/"+
							ensemblContig.apply(vc) +
							"-"+vc.getStart()+
							" "+
							vc.getReference().getDisplayString()+
							">"+
							alt.getDisplayString()
							));
						}
					for(final Allele alt: vc.getAlternateAlleles())
						{
						if(vc.getReference().isSymbolic() || alt.isSymbolic()) continue;
						//gnomad
						t.addRow("Gnomad",new HyperlinkDecorator("http://gnomad.broadinstitute.org/variant/"+
							ensemblContig.apply(vc) +
							"-"+vc.getStart()+
							"-"+
							vc.getReference().getDisplayString()+
							"-"+
							alt.getDisplayString()
							));
						}
					
					t.addRow("clinvar 37",new HyperlinkDecorator("https://www.ncbi.nlm.nih.gov/clinvar/?term="+
							ensemblContig.apply(vc) +
							"%5Bchr%5D+AND+"+ vc.getStart()+"%3A"+vc.getEnd()+"%5Bchrpos37%5D"
							));
					if(vc.getStart()!=vc.getEnd()) {
						t.addRow("decipher",new HyperlinkDecorator("https://decipher.sanger.ac.uk/search?q="+
								vc.getContig() + 
								"%3A"+vc.getStart()+"-"+vc.getEnd()));
						
						t.addRow("dgv",new HyperlinkDecorator("http://dgv.tcag.ca/gb2/gbrowse/dgv2_hg19?name="+
								ucscContig.apply(vc) + 
								"%3A"+vc.getStart()+"-"+vc.getEnd() + ";search=Search"
								)
								);

						
						String build="";
						if(SequenceDictionaryUtils.isGRCh37(header)) build="hg19";
						if(SequenceDictionaryUtils.isGRCh38(header)) build="hg38";
						if(!StringUtil.isBlank(build)) {
							t.addRow("Hi-C",new HyperlinkDecorator("http://promoter.bx.psu.edu/hi-c/view.php?method=Hi-C&species=human&assembly="+build+
									"&source=inside&tissue=GM12878&type=Lieberman-raw&resolution=25&c_url=&transfer=&gene=&chr="+vc.getContig()+"&start="+vc.getStart()+"&end="+vc.getEnd()+"&sessionID=&browser=none"));	
							}
						}
					}
				this.writeTable(margin, t);
				}
			
			if(!getOwner().hideInfo)
				{		
				/* INFO */
				final Table t=new Table("key","Index","Value").setCaption("INFO");
				final Map<String,Object> atts = vc.getAttributes();
				for(final String key: new TreeSet<>(atts.keySet()))
					{
					if(key.equals(this.vcfTools.getVepPredictionParser().getTag()) && this.vcfTools.getVepPredictionParser().isValid()) continue;
					if(key.equals(this.vcfTools.getAnnPredictionParser().getTag()) && this.vcfTools.getAnnPredictionParser().isValid()) continue;
					Object v= atts.get(key);
					final List<?> L;
					if(v instanceof List)
						{
						L=(List<?>)v;
						}
					else if(v.getClass().isArray())
						{
						Object a[]=(Object[])v;
						L=Arrays.asList(a);
						}
					else
						{
						L=Collections.singletonList(v);
						}
					for(int x=0;x< L.size();++x)
						{
						t.addRow(key,(L.size()==1?null:x+1),L.get(x));
						}
					}
				this.writeTable(margin, t);
				}
					
				/** VEP */
				if(!getOwner().hidePredictions && this.vcfTools.getVepPredictionParser().isValid())
					{
					final List<String> cats = new ArrayList<>(this.vcfTools.getVepPredictionParser().getCategories());
					final Table t = new Table(cats).setCaption("VEP");
					for(VepPrediction pred: this.vcfTools.getVepPredictionParser().getPredictions(vc))
						{
						final List<Object> row =new ArrayList<>(cats.size());//signe cat's eyes nanana
						for(final String cat:cats) {
							final String valuestr = pred.get(cat);
							final Object o;
							if(StringUtil.isBlank(valuestr))
								{
								o=valuestr;
								}
							else if(cat.equals("Existing_variation") || cat.equals("RefSeq") || cat.equals("Feature") || cat.equals("Gene") || cat.equals("ENSP"))
								{
								o = new HyperlinkDecorator(valuestr);
								}
							else if(cat.equals("Consequence"))
								{
								o = new VcfToTable.SODecorator(valuestr);
								}
							else if(cat.equals("SYMBOL")) {
								o = new GenelinkDecorator(valuestr);
								}
							else if(cat.equals("HGNC_ID")) {
								o = new HgncDecorator(valuestr);
								}
							else
								{
								o=valuestr;
								}
							row.add(o);	
							}
						t.addList(row);
						}
					t.removeEmptyColumns();
					this.writeTable(margin, t);
					}

					
				/** ANN */
				if(!getOwner().hidePredictions && this.vcfTools.getAnnPredictionParser().isValid())
					{
					Table t = new Table("SO","Allele","Impact","GeneName","GeneId","FeatureType","FeatureId",
							"BioType","HGVsc","HGVsp","Rank","cDNA-pos","CDS-pos","AA-pos","Distance","Msg").setCaption("ANN");
					
					for(final AnnPrediction P: this.vcfTools.getAnnPredictionParser().getPredictions(vc)) {
						final List<Object> r=new ArrayList<>();
						r.add(new SODecorator(P.getSOTermsString()));
						r.add(P.getAllele());
						r.add(	
									!useANSIColors || P.getPutativeImpact()==null ||P.getPutativeImpact().equals(AnnPredictionParser.Impact.UNDEFINED) || P.getPutativeImpact().equals(AnnPredictionParser.Impact.LOW)? 
									P.getPutativeImpact():
									new ColoredDecorator(P.getPutativeImpact(), AnsiColor.RED)
									);
						r.add(new GenelinkDecorator(P.getGeneName()));
						r.add(new GenelinkDecorator(P.getGeneId()));
						r.add(P.getFeatureType());
						r.add(new GenelinkDecorator(P.getFeatureId()));
						r.add(P.getTranscriptBioType());
						r.add(P.getHGVSc());
						r.add(P.getHGVSp());
						r.add(P.getRank());
						r.add(P.getCDNAPos());
						r.add(P.getCDSPos());
						r.add(P.getAAPos());
						r.add(P.getDistance());
						r.add(P.getMessages());
						t.addList(r);
						}			
					t.removeEmptyColumns();
					this.writeTable(margin, t);
					}
				
				printGenotypesTypes(margin,vc);
				printCharts(margin,vc);
				
				
				if(!getOwner().hideGenotypes && vc.hasGenotypes())
					{
					//margin = margin+ DEFAULT_MARGIN;
					final CharSplitter tab = CharSplitter.TAB;
					final CharSplitter colon = CharSplitter.COLON;
					final List<String> hds = new ArrayList<>();
					
					hds.add("Sample");
					hds.add("Type");
					
					final int prefix_header_size = hds.size();
					
					hds.addAll(header.getFormatHeaderLines().
							stream().
							map(F->F.getID()).
							collect(Collectors.toList())
							);
					final Table t=new Table(hds).setCaption("Genotypes");
					final String tokens[]=tab.split(this.vcfEncoder.encode(vc));
					final List<String> formats=Arrays.asList(colon.split(tokens[8]));
					for(int i=0;i< vc.getNSamples();i++)
						{
						final Genotype g=vc.getGenotype(i);
						if(getOwner().hideHomRefGenotypes && g.isHomRef()) continue;
						if(getOwner().hideNoCallGenotypes && !g.isCalled()) continue;
						
						final List<String> gstr =Arrays.asList(colon.split(tokens[9+i]));
						final List<Object> r= new ArrayList<>(hds.size());
						r.add(g.getSampleName());
						if(!useANSIColors ||
							g.getType().equals(GenotypeType.NO_CALL) || 
							g.getType().equals(GenotypeType.UNAVAILABLE)
							)
							{
							r.add(g.getType().name());
							}
						else if(g.getType().equals(GenotypeType.HOM_REF))
							{
							r.add(new ColoredDecorator(g.getType().name(),AnsiColor.GREEN));
							}
						else if(g.getType().equals(GenotypeType.HET))
							{
							r.add(new ColoredDecorator(g.getType().name(),AnsiColor.YELLOW));
							}
						else if(g.getType().equals(GenotypeType.HOM_VAR))
							{
							r.add(new ColoredDecorator(g.getType().name(),AnsiColor.RED));
							}
						else
							{
							r.add(new ColoredDecorator(g.getType().name(),AnsiColor.MAGENTA));
							}
						
						for(int j=prefix_header_size;j< hds.size();++j)
							{
							int indexInFORMAT = formats.indexOf(hds.get(j));
							if( indexInFORMAT==-1 || indexInFORMAT>=gstr.size()) {
								r.add(null);
								}
							else if(useANSIColors && g.isCalled() && formats.get(indexInFORMAT).equals("GQ") && g.hasGQ())  {
								final int gq = g.getGQ();
								if(gq==0)
									{
									r.add(new ColoredDecorator(String.valueOf(gq),AnsiColor.RED));
									}
								else
									{
									r.add(String.valueOf(gq));
									}
								}
							else if(useANSIColors && g.isCalled() && formats.get(indexInFORMAT).equals("DP") && g.hasDP())  {
								final int dp = g.getDP();
								r.add(new ColoredDecorator(String.valueOf(dp),dp >= 10?AnsiColor.GREEN: AnsiColor.RED));
								}
							else if(useANSIColors && g.isCalled() && formats.get(indexInFORMAT).equals("AD") && g.hasAD())  {
								final int ad_array[] = g.getAD();
								final AnsiColor ad_color;
								
								if(ad_array.length==2 && g.isHet()) {
									double ratio =  ad_array[0]/(float)(ad_array[0]+ad_array[1]);
									ad_color = ratio>=0.3 && ratio <=0.7 ? AnsiColor.GREEN: AnsiColor.RED;
									}
								else if(ad_array.length==2 && g.isHomRef())
									{
									ad_color = ad_array[1]==0 && ad_array[0]>0 ? AnsiColor.GREEN: AnsiColor.RED;
									}
								else if(ad_array.length==2 && g.isHomVar())
									{
									ad_color = ad_array[1]>0 && ad_array[0]==0 ? AnsiColor.GREEN: AnsiColor.RED;
									}
								else
									{
									ad_color = null;
									}
									
								
								if(ad_color!=null) {
									r.add(new ColoredDecorator(gstr.get(indexInFORMAT),ad_color));
									}
								else
									{
									r.add(gstr.get(indexInFORMAT));
									}
								}
							else
								{
								r.add(gstr.get(indexInFORMAT));
								}
							}
						t.addList(r);
						}
					t.removeEmptyColumns();
					this.writeTable(margin, t);
					
					
					
					if(this.pedigree!=null) {
						final Function<Genotype,String> genotype2str = G->
							G.getAlleles().stream().
								map(A->A.isNoCall()?Allele.NO_CALL_STRING:String.valueOf(vc.getAlleleIndex(A))).
								collect(Collectors.joining(G.isPhased()?"|":"/"))
							;
							
						final Table t2=new Table(
								"Father-ID",
								"Father-GT",
								"Mother-ID",
								"Mother-GT",
								"Child-ID",
								"Child-GT",
								"Incompat."
								).setCaption("TRIOS");
						for(final String childId:this.header.getSampleNamesInOrder())
							{
							final Pedigree.Person child = this.pedigree.getPersonById(childId);
							if(child==null) continue;
							final Genotype gc = vc.getGenotype(childId);
							if(gc==null) continue;
							
							final  Pedigree.Person father= child.getFather();
							final Genotype gf =  (father==null?null:vc.getGenotype(father.getId()));
							
							final  Pedigree.Person mother= child.getMother();
							final Genotype gm =  (mother==null?null:vc.getGenotype(mother.getId()));
							
							if(gf==null && gm==null) continue;
							
							final boolean is_incompat= this.vcfTools.isMendelianIncompatibility(gc, gf, gm);
							final Function<Object, Object> colorize = O->{
								if(O==null || !is_incompat|| !useANSIColors) return O;
								return new ColoredDecorator(O, AnsiColor.RED);
								};
							
							final List<Object> r= new ArrayList<>();
							r.add(father==null?null:colorize.apply(father.getId()));
							r.add(gf==null?null:colorize.apply(genotype2str.apply(gf)));
							r.add(mother==null?null:colorize.apply(mother.getId()));
							r.add(gm==null?null:colorize.apply(genotype2str.apply(gm)));
							r.add(colorize.apply(child.getId()));
							r.add(colorize.apply(genotype2str.apply(gc)));
							r.add(is_incompat?colorize.apply("*"):null);
							
							
							t2.addList(r);
							}
						this.writeTable(margin, t2);
						
					}
					
					
					
					}
				this.endVariant(vc);
				
						
				}
			
			protected String variantToString(final VariantContext vc)
				{
				final Optional<String> buildLabel = SequenceDictionaryUtils.getBuildName( header.getSequenceDictionary());
				
				return (buildLabel.isPresent()?buildLabel.get()+" ":"")+ vc.getContig()+":"+vc.getStart()+(vc.getStart()!=vc.getEnd()?"-"+vc.getEnd():"")+"/"+vc.getReference().getDisplayString();
				}
			
			protected void printGenotypesTypes(final String margin,final VariantContext vc) {
				if(getOwner().hideGTypes) return;
				if(!vc.hasGenotypes()) return;
				final Table t=new Table("Type","Count","%").
						setCaption("Genotype Types");
				vc.getGenotypes().stream().map(G->G.getType()).
					collect(Collectors.groupingBy(  Function.identity(), Collectors.counting())).
						entrySet().
						stream().
						sorted((A,B)->B.getValue().compareTo(A.getValue())).
						map(E->new Object[] {(Object)E.getKey().name(),(Object)E.getValue(),(Object)Integer.valueOf((int)(100.0*(E.getValue()/(double)vc.getNSamples())))}).
						forEach(R->t.addRow(R))
						;
				
				t.removeEmptyColumns();
				this.writeTable(margin, t);
				}
			protected void printCharts(final String margin,final VariantContext vc) {
				
				}
			}
		
		private class TerminalViewer extends AbstractViewer
			{
			private PrintStream out= VcfToTableViewer.this.outputStream;
			
			@Override
			void writeTable(final String margin,final Table t) {
				t.print(margin,this.out);
				}
			
			@Override
			void println(String s)
				{
				this.out.println(s);
				}
			
			@Override
			void println() {
				out.println();
				}

			
			
			
			@Override
			void startVariant(final VariantContext ctx) {
				this.println(">>"+ variantToString(ctx) +" (n. "+countVariants+")");
				}
			@Override
			void endVariant(final VariantContext ctx) {
				this.println("<<"+variantToString(ctx)+" (n. "+countVariants+")");				
				}
			
			@Override
			public void setHeader(final VCFHeader header) {
				throw new IllegalStateException("setHeader shouldn't be called"); 
				}
			
			@Override
			public void writeHeader(final VCFHeader header)
				{
				if(getOwner().outputFile!=null) {
					try {
						this.out = new PrintStream(IOUtils.openFileForWriting(getOwner().outputFile));
					} catch (final IOException e) {
						throw new RuntimeIOException(e);
						}
					
					}
				super.writeHeader(header);
				}
			
			
			@Override
			public boolean checkError() {
				if(out==null) return true;
				return out.checkError();
				}
			@Override
			public void close() {
				if(out==null) return;
				out.flush();
				out.close();
				out=null;
				}
			}

		private class HtmlViewer extends AbstractViewer
			{
			XMLStreamWriter out=null;
			
			protected String getCssStyle() {
				return DEFAULT_CSS_STYLE;
				}
			
			@Override
			public void setHeader(final VCFHeader header) {
				throw new IllegalStateException("setHeader shouldn't be called"); 
				}
			
			@Override
			public void writeHeader(final VCFHeader header)
				{
				try {
					final XMLOutputFactory xof=XMLOutputFactory.newFactory();
					if(getOwner().outputFile==null)
						{
						this.out = xof.createXMLStreamWriter(VcfToTableViewer.this.outputStream, "UTF-8");
						}
					else
						{
						this.out = xof.createXMLStreamWriter(new StreamResult(getOwner().outputFile));
						}
					
					if(!getOwner().hideHtmlHeader) {
						out.writeStartElement("html");
						out.writeStartElement("head");
						
						out.writeStartElement("title");
						out.writeCharacters("VcfToTable");
						out.writeEndElement();//title
						
						out.writeStartElement("style");
						out.writeCharacters(getCssStyle());
						out.writeEndElement();//style
						
						if(getOwner().googleChart)
							{
							out.writeStartElement("script");
							out.writeAttribute("type", "text/javascript");
							out.writeAttribute("src", "https://www.gstatic.com/charts/loader.js");
							out.writeCharacters("");//force blank
							out.writeEndElement();//script
							//
							out.writeStartElement("script");
							out.writeAttribute("type", "text/javascript");
							out.writeCharacters("google.charts.load('current', {'packages':['corechart']});");
							out.writeEndElement();//script
							}
						
						out.writeEndElement();//head
						out.writeStartElement("body");
						}
					this.out.writeComment("BEGIN-VCF");
					this.out.writeStartElement("div");//main-div
					this.out.writeStartElement("div");//header
					this.out.writeComment("BEGIN VCF HEADER");
					}
				catch(XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}
				
				super.writeHeader(header);
				try {
					this.out.writeComment("END VCF HEADER");
					this.out.writeEndElement();//header
					this.out.writeStartElement("div");//div-variants
					}
				catch(XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}
				}
			@Override
			public void close()
				{
				if(this.out==null) return;
				try {
					
					
					
					this.out.writeEndElement();//div-variants
					this.out.writeEndElement();//main-div
					
					this.out.writeComment("END-TABLE");
					this.out.writeComment("END-VCF");
					if(!getOwner().hideHtmlHeader)
						{
						this.out.writeEndElement();//body
						this.out.writeEndElement();//html
						}
					
					this.out.flush();
					this.out.close();
					}
				catch(final XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}
				out=null;
				}
			@Override
			void println()
				{
				try {
					this.out.writeEmptyElement("br");
					}
				catch(final XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}
				}
			@Override
			void println(final String s)
				{
				try {
					this.out.writeStartElement("p");
					this.out.writeCharacters(s==null?"":s);
					this.out.writeEndElement();
					}
				catch(final XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}
				}
			@Override
			void writeTable(String margin, final Table t)
				{
				try {
					
					t.write(this.out);
					}
				catch(final XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}
				}
			@Override
			public boolean checkError()
				{
				if(this.out==null) return true;
				return false;
				}
			
			@Override
			void startVariant(final VariantContext ctx) {
				try {
					this.out.writeStartElement("div");//div for whole variant
					this.out.writeAttribute("class","variant"+(this.countVariants%2));
					
					this.out.writeEmptyElement("a");
					this.out.writeAttribute("name","vc"+this.countVariants);
					this.out.writeStartElement("h3");
					this.out.writeCharacters(variantToString(ctx)+" (n. "+this.countVariants+").");
					if(this.countVariants>1)
						{
						this.out.writeStartElement("a");
						this.out.writeAttribute("href", "#vc"+(this.countVariants-1));
						this.out.writeCharacters("[prev]");
						this.out.writeEndElement();//a
						}
					
					/* Hyperlink to IGV */
					if(!StringUtil.isBlank(VcfToTableViewer.this.userCustomUrl)) {
						final int EXTEND=15;
						final String gotostr=Launcher.createUrlFromInterval(
								VcfToTableViewer.this.userCustomUrl,
								new Interval(ctx.getContig(),
									Math.max(1, ctx.getStart()-EXTEND),
									ctx.getEnd()+EXTEND));
						if(!StringUtil.isBlank(gotostr)) {
							this.out.writeCharacters(" ");
							this.out.writeStartElement("a");
							this.out.writeAttribute("title","Open URL");
							this.out.writeAttribute("rel","nofollow");
							this.out.writeAttribute("href", gotostr );
							this.out.writeCharacters("[URL]");
							this.out.writeEndElement();//a
							}
						}
					
					this.out.writeCharacters(" ");
					this.out.writeStartElement("a");
					this.out.writeAttribute("href", "#vc"+(this.countVariants+1));
					this.out.writeCharacters("[next]");
					this.out.writeEndElement();//a
					
					
					
					this.out.writeEndElement();//h3
					
					
					
					}
				catch(final XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}
				
				}
			@Override
			void endVariant(final VariantContext ctx) {
				try {
					out.writeEmptyElement("br");
					this.out.writeCharacters(variantToString(ctx)+" (n. "+this.countVariants+"). ");
					
					this.out.writeStartElement("a");
					this.out.writeAttribute("href", "#vc"+(this.countVariants));
					this.out.writeCharacters("[top]");
					this.out.writeEndElement();//a
					
					
					this.out.writeEndElement();//div for variant
					this.out.writeEmptyElement("hr");
					this.out.flush();
					}
				catch(final XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}				
				}
			@Override
			protected void printGenotypesTypes(String margin, VariantContext vc) {
				if(!getOwner().googleChart) {
					super.printGenotypesTypes(margin, vc);
					return;
					}
				if(getOwner().hideGTypes) return;
				if(!vc.hasGenotypes()) return;
				//  <div id="piechart" style="width: 900px; height: 500px;"></div>
				try {
					final String id = "gtypes"+countVariants;
					final StringWriter strw = new StringWriter();
					final PrintWriter pw= new PrintWriter(strw);
					pw.print("google.charts.setOnLoadCallback(function(){");
					pw.print("var data = google.visualization.arrayToDataTable([");
					pw.print("['Type', 'Count']");
					for(final GenotypeType gt:GenotypeType.values())
						{
						final long n = vc.getGenotypes().stream().filter(G->G.getType().equals(gt)).count();
						if(n==0L) continue;
						pw.print(",['"+gt.name()+"',"+n+"]");
						}
					pw.print("]);");
					pw.print("var chart = new google.visualization.PieChart(document.getElementById('"+id+"'));");
					pw.print("chart.draw(data, {\"title\":\"Genotype Types\"});");
					pw.print("});");
					pw.flush();
					
					
					this.out.writeStartElement("span");
					this.out.writeAttribute("id", id);
					this.out.writeEndElement();
					this.out.writeStartElement("script");
					this.out.writeCharacters(strw.toString());
					this.out.writeEndElement();
					}
				catch(final XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}	
				}
			
			protected void printCharts(final String margin,final VariantContext vc) {
				if(!getOwner().googleChart) return;
				try {
					printIntChart(vc,"dp","DP",V->V.hasDP(),V->V.getDP());
					printIntChart(vc,"gq","GQ",V->V.hasGQ(),V->V.getGQ());
					printADChart(vc,"adx","HOM_REF",G->G.isHomRef() && G.hasAD() && Arrays.stream(G.getAD()).skip(1L).anyMatch(I->I>0));
					printADChart(vc,"ady","HET",G->G.isHet() && G.hasAD() );
					printADChart(vc,"adz","HOM_VAR",G->G.isHomVar() && G.hasAD() && Arrays.stream(G.getAD()).limit(1L).anyMatch(I->I>0));
					}
				catch(final XMLStreamException err)
					{
					throw new RuntimeIOException(err);
					}
				}
			private void printIntChart(
						final VariantContext vc,
						final String prefix,
						final String title,
						final Predicate<Genotype> hasInfo,
						final Function<Genotype,Integer> extractor
						) throws XMLStreamException {
				if(!vc.hasGenotypes()) return;
				if(vc.getGenotypes().stream().noneMatch(hasInfo)) return;
				
				final String id = prefix + countVariants;
				final StringWriter strw = new StringWriter();
				final PrintWriter pw= new PrintWriter(strw);
				pw.print("google.charts.setOnLoadCallback(function(){");
				pw.print("var data = google.visualization.arrayToDataTable([");
				pw.print("[\"Sample\",\""+title+"\"]");
				for(final Genotype gt : vc.getGenotypes().
						stream().
						filter(hasInfo).
						sorted((A,B)->extractor.apply(A).compareTo(extractor.apply(B))).
						collect(Collectors.toList()))
					{
					pw.print(",['"+gt.getSampleName()+"',"+extractor.apply(gt)+"]");
					}
				pw.print("]);");
				pw.print("var chart = new google.visualization.ColumnChart(document.getElementById('"+id+"'));");
				pw.print("chart.draw(data, {"+getGoogleChartDimension()+",\"title\":\""+title+"\",\"legend\": {\"position\":\"none\" }});");
				pw.print("});");
				pw.flush();
				this.out.writeStartElement("span");
				this.out.writeAttribute("id", id);
				this.out.writeEndElement();
				this.out.writeStartElement("script");
				this.out.writeCharacters(strw.toString());
				this.out.writeEndElement();
				}
			
				private void printADChart(
						final VariantContext vc,
						final String prefix,
						final String title,
						final Predicate<Genotype> hasType
						) throws XMLStreamException {
				if(!vc.hasGenotypes()) return;
				if(!vc.isVariant()) return;
				if(vc.getGenotypes().stream().filter(G->G.hasAD()).filter(hasType).noneMatch(V->V.hasAD())) return;
				
				final String id = prefix+countVariants;
				final StringWriter strw = new StringWriter();
				final PrintWriter pw= new PrintWriter(strw);
				pw.print("google.charts.setOnLoadCallback(function(){");
				pw.print("var data = google.visualization.arrayToDataTable([");
				pw.print("[\"Sample\"");
				for(int x=0;x< vc.getNAlleles();++x) pw.print(",\"Allele "+x+"\"");
				pw.print("]");
				for(final Genotype gt:vc.getGenotypes().
						stream().
						filter(G->G.hasAD()).
						filter(hasType).
						sorted((A,B)->Integer.compare(A.getAD()[0], B.getAD()[0])).
						collect(Collectors.toList()))
					{
					final int array[] = gt.getAD();
					final double sum = IntStream.of(array).sum();
					if(sum==0.0) continue;
					pw.print(",['"+gt.getSampleName()+"'");
					for(int x=0;x< vc.getNAlleles();++x)
						{
						pw.print(",");
						pw.print((int)(((x<array.length?array[x]:0)/sum)*100.0));
						}
					pw.print("]");
					}
				pw.print("]);");
				pw.print("var chart = new google.visualization.ColumnChart(document.getElementById('"+id+"'));");
				pw.print("chart.draw(data, {"+getGoogleChartDimension()+",\"isStacked\": true,\"title\":\"AD per sample ("+title+")\",\"legend\": {\"position\":\"bottom\" }});");
				pw.print("});");
				pw.flush();
				this.out.writeStartElement("span");
				this.out.writeAttribute("id", id);
				this.out.writeEndElement();
				this.out.writeStartElement("script");
				this.out.writeCharacters(strw.toString());
				this.out.writeEndElement();
				}
			}

		
		/** public, to be used in e.g; a VCF server */
		public VcfToTableViewer() {
			}
		
		/** give a chance to the class using vcf2table to change output stream */
		public void setOutputStream(final PrintStream outputStream) {
			this.outputStream = outputStream;
			}
		
		public void setOutputFormat(final OutputFormat outputFormat) {
			this.outputFormat = outputFormat;
			}
		
		public void setHideHtmlHeader(boolean hideHtmlHeader) {
			this.hideHtmlHeader = hideHtmlHeader;
			}
		
		public void setPrintHeader(boolean printHeader) {
			this.printHeader = printHeader;
		}
		
		public void setPedigreeFile(File pedigreeFile) {
			this.pedigreeFile = pedigreeFile;
			}
		
		public void setHideHomRefGenotypes(boolean hideHomRefGenotypes) {
			this.hideHomRefGenotypes = hideHomRefGenotypes;
		}
		
		public void setHideNoCallGenotypes(boolean hideNoCallGenotypes) {
			this.hideNoCallGenotypes = hideNoCallGenotypes;
		}
		public void setHideGenotypes(boolean hideGenotypes) {
			this.hideGenotypes = hideGenotypes;
		}
		public void setUseANSIColors(boolean useANSIColors) {
			this.useANSIColors = useANSIColors;
		}
		public void setUserCustomUrl(final String userCustomUrl) {
			this.userCustomUrl = userCustomUrl;
		}
		public void setHideGTypes(boolean hideGTypes) {
			this.hideGTypes = hideGTypes;
		}
		
		@Override
		public void setHeader(final VCFHeader header) {
			throw new IllegalStateException("setHeader shouldn't be called"); 
			}
		
		
		@Override
		public void writeHeader(final VCFHeader header) {
			if(this.delegate!=null) throw new IllegalStateException(); 
			switch(this.outputFormat)
				{
				case html: this.delegate = new HtmlViewer();break;
				default: this.delegate = new TerminalViewer();break;
				}
			
			this.delegate.writeHeader(header);
			}
		@Override
		public void add(final VariantContext vc) {
			if(this.delegate==null) return ; 
			this.delegate.add(vc);
			}
		@Override
		public boolean checkError() {
			if(this.delegate==null) return true;
			return this.delegate.checkError();
			}
		@Override
		public void close() {
			if(this.delegate==null) return;
			this.delegate.close();
			this.delegate=null;
			}
		
		
		String getGoogleChartDimension() {
			int w=0,h=0;
			final int x= StringUtil.isBlank(this.googleChartSizeStr)?-1:this.googleChartSizeStr.indexOf('x');
			try {
				if(x>1) {
					w = Integer.parseInt(this.googleChartSizeStr.substring(0, x));
					h = Integer.parseInt(this.googleChartSizeStr.substring(x+1));
					}
				else if(!StringUtil.isBlank(this.googleChartSizeStr))
					{
					w = Integer.parseInt(this.googleChartSizeStr);
					h = (int)(w/1.618);
					}
				}
			catch(final NumberFormatException err) {
				w = -1;
				h = -1;
				}
			if(!(w>0 && h>0)) { w=1000;h=618;}
			return  "\"width\":"+w+",\"height\":"+h;
			}

		}
	
	
	@ParametersDelegate
	private VcfToTableViewer  viewer = new VcfToTableViewer();

	
	
	@Override
	public int doWork(final List<String> args) {
		VCFIterator in = null;
		
		
		try {
			in = super.openVCFIterator(oneFileOrNull(args));
			viewer.writeHeader(in.getHeader());
			while(!this.viewer.checkError() && in.hasNext())
				{
				viewer.add(in.next());
				}
			viewer.close();viewer=null;
			in.close();in=null;
			return 0;
		} catch (Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(viewer);
			}
		}
	public static void main(final String[] args) {
		new VcfToTable().instanceMainWithExit(args);
	}
}
