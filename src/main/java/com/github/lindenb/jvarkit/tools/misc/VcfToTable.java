/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
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
END_DOC

## History

* 20170517 : fixed a bug in removeEmptyColumns

 */
@Program(name="vcf2table",
		description="convert a vcf to a table, to ease display in the terminal",
		keywords={"vcf","table","visualization"})
public class VcfToTable extends Launcher {
	private static final Logger LOG = Logger.build(VcfToTable.class).make();
	private static final String DEFAULT_MARGIN=" ";
	public static final String ANSI_ESCAPE = "\u001B[";
	public static final String ANSI_RESET = ANSI_ESCAPE+"0m";
	
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

	private static class Colored
		{
		final Object o;
		final AnsiColor color;
		Colored(final Object o,final AnsiColor color) { this.o=o;this.color=color;}
		public String start() {
			if(o==null || color==null) return "";
			return ANSI_ESCAPE+this.color.opcode+"m";
			}
		public String end()
			{
			if(o==null || color==null) return "";
			return ANSI_RESET;
			}
		
		@Override
		public String toString() {
			return o==null?"":o.toString();
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
					if(cell instanceof Colored) out.print(Colored.class.cast(cell).start());
					out.print(str);
					if(cell instanceof Colored) out.print(Colored.class.cast(cell).end());
					for(int j=str.length();j< this.columns.get(i).maxLength;j++) out.print(' ');
					out.print(" |");
					}
				out.println();
				
			}
			out.print(margin);
			out.println(hr.toString());
					
			}
		}
	
	public static class TerminalViewer
		implements VariantContextWriter
		{
		@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
		private File outputFile = null;
		@Parameter(names={"-H"},description="Print Header")
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
		
		private int countVariants=0;
		
		private PrintStream out= System.out;
		private VCFHeader header=null;
		private VCFEncoder vcfEncoder=null;
		private Pedigree pedigree = null;
		private VcfTools vcfTools = null;
		TerminalViewer() {
			}
		@Override
		public void writeHeader(final VCFHeader header) {
			String margin="";
			this.header = header;
			this.vcfEncoder = new VCFEncoder(header, true, true);
			this.vcfTools = new VcfTools(header);
			if(outputFile!=null) {
				try {
					this.out = new PrintStream(IOUtils.openFileForWriting(this.outputFile));
				} catch (final IOException e) {
					throw new RuntimeIOException(e);
					}
				
				}
			
			if(this.pedigreeFile!=null) {
				try {
					this.pedigree = Pedigree.newParser().parse(this.pedigreeFile);
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
			
			if(printHeader)
				{
				/** INFO */
					{
					Table t=new Table("ID","Type","Count","Description").setCaption("INFO");
					header.getInfoHeaderLines().stream().
						map(F->{
								final List<Object> r=new ArrayList<>();
								r.add(F.getID());
								r.add(F.getType()==null?null:F.getType().name());
								r.add(F.isFixedCount()?F.getCount():null);
								r.add(F.getDescription());
								return r;
								}).
						forEach(R->t.addList(R));
					t.print(margin,out);
					out.println();
					}
				/** FORMAT */
					{
					Table t=new Table("ID","Type","Count","Description").setCaption("FORMAT");
					header.getFormatHeaderLines().stream().
						map(F->{
								final List<Object> r=new ArrayList<>();
								r.add(F.getID());
								r.add(F.getType()==null?null:F.getType().name());
								r.add(F.isFixedCount()?F.getCount():null);
								r.add(F.getDescription());
								return r;
								}).
						forEach(R->t.addList(R));
					t.print(margin,out);
					out.println();
					}
				/** FILTER */
					{
					final Table t=new Table("ID","Description").setCaption("FILTERS");
					header.getFilterLines().forEach(
						L->t.addRow(L.getID(),L.getDescription())
						);
					
					t.removeEmptyColumns().print(margin,out);
					out.println();
					}
				
				/** OTHER METADATA */
					{
					final Table t=new Table("ID","Description").setCaption("Metadata");
					header.getOtherHeaderLines().forEach(
						L->t.addRow(L.getKey(),L.getValue())
						);
					
					t.removeEmptyColumns().print(margin,out);
					out.println();
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
						t2.print(margin,out);
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
					
					t.print(margin,out);
					}
					
					
				out.println();
				}
			}
		@Override
		public void add(final VariantContext vc) {
			if(out==null) return;
			
			if(this.limitVariants!=-1)
				{
				if(this.countVariants==this.limitVariants)
					{
					out.println();
					out.println("... More variants exist but they've been omitted...( limit reached)");
					return;
					}
				else if(this.countVariants>this.limitVariants)
					{
					return;
					}
				}
			
			final String variantName=vc.getContig()+"/"+vc.getStart()+"/"+vc.getReference().getDisplayString();
			++countVariants;
			out.println(">>"+ variantName+" (n. "+countVariants+")");
			String margin=DEFAULT_MARGIN;
			{
			final Table t=new Table("Key","Value").setCaption("Variant");
			t.addRow("CHROM",vc.getContig());
			t.addRow("POS",vc.getStart());
			t.addRow("end",vc.getEnd());
			t.addRow("ID",vc.hasID()?vc.getID():".");
			t.addRow("REF",vc.getReference().getDisplayString());
			t.addRow("ALT",vc.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(",")));
			t.addRow("QUAL",vc.hasLog10PError()?vc.getPhredScaledQual():null);
			final String filterStr=vc.isFiltered()?vc.getFilters().stream().collect(Collectors.joining(";")):null;
			t.addRow("FILTER",filterStr!=null && useANSIColors?new Colored(filterStr, AnsiColor.RED):filterStr);
			t.addRow("Type",vc.getType());
			
			
			
			t.print(margin,out);
			}
		
		if(!this.hideAlleles)
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
			t.print(margin,out);
			}
		
		if(!this.hideFilters)
			{
			/* FILTER */
			final	Table t=new Table("Filter").setCaption("FILTERS");
			 for(final String f:vc.getFilters())
			 	{
				t.addRow(this.useANSIColors?new Colored(f, AnsiColor.YELLOW):f);
			 	}
			t.print(margin,out);
			}
		
		if(!this.hideInfo)
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
			t.print(margin,out);
			}
				
			/** VEP */
			if(!this.hidePredictions && this.vcfTools.getVepPredictionParser().isValid())
				{
				final List<String> cats = new ArrayList<>(this.vcfTools.getVepPredictionParser().getCategories());
				
				final Table t = new Table(cats).setCaption("VEP");
				for(VepPrediction pred: this.vcfTools.getVepPredictionParser().getPredictions(vc))
					{
					final List<Object> r=new ArrayList<>();
					for(final String cat:cats) {
						r.add(pred.get(cat));
						}
					t.addList(r);
					}
				t.removeEmptyColumns();
				t.print(margin,out);
				}

				
			/** ANN */
			if(!this.hidePredictions && this.vcfTools.getAnnPredictionParser().isValid())
				{
				Table t = new Table("SO","Allele","Impact","GeneName","GeneId","FeatureType","FeatureId",
						"BioType","HGVsc","Rank","cDNA-pos","CDS-pos","AA-pos","Distance","Msg").setCaption("ANN");
				
				for(final AnnPrediction P: this.vcfTools.getAnnPredictionParser().getPredictions(vc)) {
					final List<Object> r=new ArrayList<>();
					r.add(P.getSOTermsString());
					r.add(P.getAllele());
					r.add(	
								!useANSIColors || P.getPutativeImpact()==null ||P.getPutativeImpact().equals(AnnPredictionParser.Impact.UNDEFINED) || P.getPutativeImpact().equals(AnnPredictionParser.Impact.LOW)? 
								P.getPutativeImpact():
								new Colored(P.getPutativeImpact(), AnsiColor.RED)
								);
					r.add(P.getGeneName());
					r.add(P.getGeneId());
					r.add(P.getFeatureType());
					r.add(P.getFeatureId());
					r.add(P.getTranscriptBioType());
					r.add(P.getHGVSc());
					r.add(P.getRank());
					r.add(P.getCDNAPos());
					r.add(P.getCDSPos());
					r.add(P.getAAPos());
					r.add(P.getDistance());
					r.add(P.getMessages());
					t.addList(r);
					}			
				t.removeEmptyColumns();
				t.print(margin,out);
				}
			if(!this.hideGenotypes)
				{
				//margin = margin+ DEFAULT_MARGIN;
				final Pattern tab = Pattern.compile("\t");
				final Pattern colon = Pattern.compile("\\:");
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
					if(this.hideHomRefGenotypes && g.isHomRef()) continue;
					if(this.hideNoCallGenotypes && !g.isCalled()) continue;
					
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
						r.add(new Colored(g.getType().name(),AnsiColor.GREEN));
						}
					else if(g.getType().equals(GenotypeType.HET))
						{
						r.add(new Colored(g.getType().name(),AnsiColor.YELLOW));
						}
					else if(g.getType().equals(GenotypeType.HOM_VAR))
						{
						r.add(new Colored(g.getType().name(),AnsiColor.RED));
						}
					else
						{
						r.add(new Colored(g.getType().name(),AnsiColor.MAGENTA));
						}
					
					for(int j=prefix_header_size;j< hds.size();++j)
						{
						int indexInFORMAT = formats.indexOf(hds.get(j));
						if( indexInFORMAT==-1 || indexInFORMAT>=gstr.size()) {
							r.add(null);
							}
						else
							{
							r.add(gstr.get(indexInFORMAT));
							}
						}
					t.addList(r);
					}
				t.removeEmptyColumns();
				t.print(margin,out);
				
				
				
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
							return new Colored(O, AnsiColor.RED);
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
					
					t2.print(margin,out);
				}
				
				
				
				}
			out.println("<<"+variantName+" (n. "+countVariants+")");
			out.println();
			out.println();
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
	
	
	@ParametersDelegate
	private TerminalViewer  viewer = new TerminalViewer();

	@Override
	public int doWork(List<String> args) {
		VcfIterator in = null;
		
		
		try {
			in = super.openVcfIterator(oneFileOrNull(args));
			viewer.writeHeader(in.getHeader());
			while(!this.viewer.checkError() && in.hasNext())
				{
				viewer.add(in.next());
				}
			viewer.out.close();viewer.out=null;
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
	public static void main(String[] args) {
		new VcfToTable().instanceMainWithExit(args);
	}
}
