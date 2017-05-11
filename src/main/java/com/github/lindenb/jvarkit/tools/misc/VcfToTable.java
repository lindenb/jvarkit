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
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;

/**

BEGIN_DOC

## Example 

```
$ curl -s "https://raw.githubusercontent.com/arq5x/gemini/master/test/test.region.vep.vcf" | java -jar dist/vcf2table.jar -H

INFO
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
| ID              | Type    | Count | Description                                                                                                                                             |
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+
| AC              | Integer |       | Allele count in genotypes, for each ALT allele, in the same order as listed                                                                             |
| AF              | Float   |       | Allele Frequency, for each ALT allele, in the same order as listed                                                                                      |
| AN              | Integer | 1     | Total number of alleles in called genotypes                                                                                                             |
| BaseQRankSum    | Float   | 1     | Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities                                                                                       |
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
| CSQ             | String  |       | Consequence type as predicted by VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|ALLELE_NUM |
+-----------------+---------+-------+---------------------------------------------------------------------------------------------------------------------------------------------------------+

INFO
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
| chr10                 | 135534747 | hg19 |
| chrY                  | 59373566  | hg19 |
+-----------------------+-----------+------+

>>chr1/10001/T 1
Variant
+-------+--------------------+
| Key   | Value              |
+-------+--------------------+
| CHROM | chr1               |
| POS   | 10001              |
| end   | 10001              |
| ID    | .                  |
| REF   | T                  |
| ALT   | TC                 |
| QUAL  | 175.91000000000003 |
+-------+--------------------+
Alleles
+-----+-----+-----+-------+--------+
| Idx | REF | Sym | Bases | Length |
+-----+-----+-----+-------+--------+
| 0   | *   |     | T     | 1      |
| 1   |     |     | TC    | 2      |
+-----+-----+-----+-------+--------+
INFO
+----------------+-------+----------------------------------------------------------------------------------------------------------+
| key            | Index | Value                                                                                                    |
+----------------+-------+----------------------------------------------------------------------------------------------------------+
| AC             |       | 4                                                                                                        |
| AF             |       | 0.50                                                                                                     |
| AN             |       | 8                                                                                                        |
| BaseQRankSum   |       | 4.975                                                                                                    |
| CSQ            | 1     | upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000456328|||||processed_transcript|1               |
| CSQ            | 2     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000488147|||||unprocessed_pseudogene|1            |
| CSQ            | 3     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000541675|||||unprocessed_pseudogene|1            |
| CSQ            | 4     | upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000450305|||||transcribed_unprocessed_pseudogene|1 |
| CSQ            | 5     | upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000515242|||||transcribed_unprocessed_pseudogene|1 |
| CSQ            | 6     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000538476|||||unprocessed_pseudogene|1            |
| CSQ            | 7     | upstream_gene_variant|||ENSG00000223972|DDX11L1|ENST00000518655|||||transcribed_unprocessed_pseudogene|1 |
| CSQ            | 8     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000438504|||||unprocessed_pseudogene|1            |
| CSQ            | 9     | downstream_gene_variant|||ENSG00000227232|WASH7P|ENST00000423562|||||unprocessed_pseudogene|1            |
| DP             |       | 76                                                                                                       |
| FS             |       | 12.516                                                                                                   |
| HRun           |       | 0                                                                                                        |
| HaplotypeScore |       | 218.6157                                                                                                 |
| MQ             |       | 35.31                                                                                                    |
| MQ0            |       | 0                                                                                                        |
| MQRankSum      |       | -0.238                                                                                                   |
| QD             |       | 2.31                                                                                                     |
| ReadPosRankSum |       | 2.910                                                                                                    |
+----------------+-------+----------------------------------------------------------------------------------------------------------+
VEP
+----------+------+------+------------+-----------------+---------+------------------+-------------------------+-------------+--------+-----------------+------------------------------------+
| PolyPhen | EXON | SIFT | ALLELE_NUM | Gene            | SYMBOL  | Protein_position | Consequence             | Amino_acids | Codons | Feature         | BIOTYPE                            |
+----------+------+------+------------+-----------------+---------+------------------+-------------------------+-------------+--------+-----------------+------------------------------------+
|          |      |      | 1          | ENSG00000223972 | DDX11L1 |                  | upstream_gene_variant   |             |        | ENST00000456328 | processed_transcript               |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000488147 | unprocessed_pseudogene             |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000541675 | unprocessed_pseudogene             |
|          |      |      | 1          | ENSG00000223972 | DDX11L1 |                  | upstream_gene_variant   |             |        | ENST00000450305 | transcribed_unprocessed_pseudogene |
|          |      |      | 1          | ENSG00000223972 | DDX11L1 |                  | upstream_gene_variant   |             |        | ENST00000515242 | transcribed_unprocessed_pseudogene |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000538476 | unprocessed_pseudogene             |
|          |      |      | 1          | ENSG00000223972 | DDX11L1 |                  | upstream_gene_variant   |             |        | ENST00000518655 | transcribed_unprocessed_pseudogene |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000438504 | unprocessed_pseudogene             |
|          |      |      | 1          | ENSG00000227232 | WASH7P  |                  | downstream_gene_variant |             |        | ENST00000423562 | unprocessed_pseudogene             |
+----------+------+------+------------+-----------------+---------+------------------+-------------------------+-------------+--------+-----------------+------------------------------------+
Genotypes
+---------+------+-------+----+-------+-----+---------+
| Sample  | Type | AD    | DP | GQ    | GT  | PL      |
+---------+------+-------+----+-------+-----+---------+
| M10475  | HET  | 10,2  | 15 | 10.41 | 0/1 | 25,0,10 |
| M10478  | HET  | 10,4  | 16 | 5.30  | 0/1 | 40,0,5  |
| M10500  | HET  | 10,10 | 21 | 7.48  | 0/1 | 111,0,7 |
| M128215 | HET  | 15,5  | 24 | 0.26  | 0/1 | 49,0,0  |
+---------+------+-------+----+-------+-----+---------+
<<1
(...)
```
END_DOC

 */
@Program(name="vcf2table",
		description="convert a vcf to a table, to ease display in the terminal",
		keywords={"vcf","table","visualization"})
public class VcfToTable extends Launcher {
	private static final Logger LOG = Logger.build(VcfToTable.class).make();

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
			this.columns=labels.stream().
					map(L->new Column(L)).
					collect(Collectors.toList());
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
			rows.add(row);
			for(int i=0;i< row.size()&& i< this.columns.size() ;++i)
				{
				final Object o = row.get(i);
				if(o==null) continue;
				final String str=o.toString();
				final int len = str.length();
				this.columns.get(i).maxLength = Math.max(this.columns.get(i).maxLength,len);
				}
			}
		/*
		private int tableWidth() {
			return this.columns.stream().mapToInt(C->C.maxLength +2 ).sum() 
					+ 1
					+ this.columns.size();
			}*/
	
		
		public void print(final PrintStream out) {
			if(this.rows.isEmpty()) return;
			final StringBuilder hr= new StringBuilder();
			
			out.println(this.caption);
			hr.append('+');
			for(int i=0;i< this.columns.size();i++)
				{
				hr.append("-");
				for(int j=0;j< this.columns.get(i).maxLength;j++)hr.append('-');
				hr.append("-+");
				}
			out.println(hr.toString());
			
			out.print('|');
			for(int i=0;i< this.columns.size();i++)
				{
				out.print(" ");
				out.print(this.columns.get(i).label);
				for(int j=this.columns.get(i).label.length();j< this.columns.get(i).maxLength;j++) out.print(' ');
				out.print(" |");
				}
			out.println();
			out.println(hr.toString());
			
			for(int y=0;y< this.rows.size();++y) {
				final List<Object> row= this.rows.get(y);
				out.print("|");
				for(int i=0;i< this.columns.size() && i< row.size();i++)
					{
					String str= row.get(i)==null?"":row.get(i).toString();
					out.print(" ");
					out.print(str);
					for(int j=str.length();j< this.columns.get(i).maxLength;j++) out.print(' ');
					out.print(" |");
					}
				out.println();
				
			}
			out.println(hr.toString());
					
			}
		}
	
	public static class TerminalViewer
		implements VariantContextWriter
		{
		@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
		private File outputFile = null;
		@Parameter(names={"-H"},description="Print Header")
		private boolean printHeader=false;
		@Parameter(names={"-g","--hideGenotypes"},description="Hide All genotypes")
		private boolean hideGenotypes=false;
		@Parameter(names={"-nc","--hideNoCalls"},description="Hide NO_CALL genotypes")
		private boolean hideNoCallGenotypes=false;
		@Parameter(names={"-hr","--hideHomRefs"},description="Hide HOM_REF genotypes")
		private boolean hideHomRefGenotypes=false;
		
		private int countVariants=0;
		
		private PrintStream out= System.out;
		private VCFHeader header=null;
		private VepPredictionParser vepPredictionParser=null;
		private AnnPredictionParser annPredictionParser=null;
		private VCFEncoder vcfEncoder=null;
		TerminalViewer() {
			}
		@Override
		public void writeHeader(VCFHeader header) {
			this.header = header;
			this.vcfEncoder = new VCFEncoder(header, true, true);
			this.vepPredictionParser = new VepPredictionParserFactory(this.header).get();
			this.annPredictionParser = new AnnPredictionParserFactory(this.header).get();
			if(outputFile!=null) {
				try {
					this.out = new PrintStream(IOUtils.openFileForWriting(this.outputFile));
				} catch (final IOException e) {
					throw new RuntimeIOException(e);
				}
				
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
					t.print(out);
					out.println();
					}
				/** FORMAT */
					{
					Table t=new Table("ID","Type","Count","Description").setCaption("INFO");
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
					t.print(out);
					out.println();
					}
				/** FILTER */
				
					{
					Table t=new Table("ID","Description").setCaption("FILTERS");
					header.getFilterLines().forEach(
						L->t.addRow(L.getID(),L.getDescription())
						);
					t.print(out);
					out.println();
					}
				
					
					
				final SAMSequenceDictionary dict = header.getSequenceDictionary();
				if (dict != null) {
					final List<String> h = new ArrayList<>();
					h.add("Name");
					h.add("Length");
					final Set<String> all_attributes = dict.getSequences().stream()
							.flatMap(S -> S.getAttributes().stream()).map(A -> A.getKey())
							.filter(S -> !(S.equals("Name") || S.equals("length") || S.equals("Length"))).collect(Collectors.toSet());
					h.addAll(all_attributes);
					Table t2 = new Table(h).setCaption("Dict");

					for (final SAMSequenceRecord ssr : dict.getSequences()) {
						final List<Object> r = new ArrayList<>();
						r.add(ssr.getSequenceName());
						r.add(ssr.getSequenceLength());
						for (final String key : all_attributes) {
							r.add(ssr.getAttribute(key));
						}
						t2.addList(r);
					}
					t2.print(out);
				}

				out.println();
				}
			}
		@Override
		public void add(final VariantContext vc) {
			if(out==null) return;
			++countVariants;
			out.println(">>"+vc.getContig()+"/"+vc.getStart()+"/"+vc.getReference().getDisplayString()+" "+countVariants);
			
			{
			final Table t=new Table("Key","Value").setCaption("Variant");
			t.addRow("CHROM",vc.getContig());
			t.addRow("POS",vc.getStart());
			t.addRow("end",vc.getEnd());
			t.addRow("ID",vc.hasID()?vc.getID():".");
			t.addRow("REF",vc.getReference().getDisplayString());
			t.addRow("ALT",vc.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(",")));
			t.addRow("QUAL",vc.hasLog10PError()?vc.getPhredScaledQual():null);
			
			
			t.print(out);
			}
		
			{
			 final Table t=new Table("Idx","REF","Sym","Bases","Length").setCaption("Alleles");
			 for(final Allele a: vc.getAlleles())
			 	{
				t.addRow(vc.getAlleleIndex(a),
						a.isReference()?"*":"",
						a.isSymbolic()?"*":"",
						a.getDisplayString(),
						a.isSymbolic()?null:a.length()
						);
			 	}
			t.print(out);
			}
			
			{
			/* FILTER */
			final	Table t=new Table("Filter").setCaption("FILTERS");
			 for(final String f:vc.getFilters())
			 	{
				t.addRow(f);
			 	}
			t.print(out);
			}
			
					{		
					/* INFO */
					final Table t=new Table("key","Index","Value").setCaption("INFO");
					final Map<String,Object> atts = vc.getAttributes();
					for(final String key: new TreeSet<>(atts.keySet()))
						{
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
					
					
					
				t.print(out);
				}
			
			/** VEP */
			if(this.vepPredictionParser.isValid())
				{
				final List<String> cats = new ArrayList<>(this.vepPredictionParser.getCategories());
				
				final Table t = new Table(cats).setCaption("VEP");
				for(VepPrediction pred: this.vepPredictionParser.getPredictions(vc))
					{
					final List<Object> r=new ArrayList<>();
					for(final String cat:cats) {
						r.add(pred.get(cat));
						}
					t.addList(r);
					}
				t.print(out);
				}

				
			/** ANN */
			if(this.annPredictionParser.isValid())
				{
				Table t = new Table("SO","Allele","Impact","GeneName","GeneId","FeatureType","FeatureId",
						"BioType","HGVsc","Rank","cDNA-pos","CDS-pos","AA-pos","Distance","Msg").setCaption("ANN");
				
				for(final AnnPrediction P: this.annPredictionParser.getPredictions(vc)) {
					final List<Object> r=new ArrayList<>();
					r.add(P.getSOTermsString());
					r.add(P.getAllele());
					r.add(P.getPutativeImpact());
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
				t.print(out);
				}
			if(!this.hideGenotypes)
				{
				final Pattern tab = Pattern.compile("\t");
				final Pattern colon = Pattern.compile("\\:");
				final List<String> hds = new ArrayList<>();
				hds.add("Sample");
				hds.add("Type");
				hds.addAll(header.getFormatHeaderLines().
						stream().
						map(F->F.getID()).
						collect(Collectors.toList()));
				final Table t=new Table(hds).setCaption("Genotypes");
				String tokens[]=tab.split(this.vcfEncoder.encode(vc));
				final List<String> formats=Arrays.asList(colon.split(tokens[8]));
				for(int i=0;i< vc.getNSamples();i++)
					{
					final Genotype g=vc.getGenotype(i);
					if(this.hideHomRefGenotypes && g.isHomRef()) continue;
					if(this.hideNoCallGenotypes && g.isNoCall()) continue;
					final List<String> gstr =Arrays.asList(colon.split(tokens[9+i]));
					final List<Object> r= new ArrayList<>(hds.size());
					r.add(g.getSampleName());
					r.add(g.getType().name());
					for(int j=2 /* 0== sample 1==type*/;j< hds.size();++j)
						{
						int x=formats.indexOf(hds.get(j));
						if( x==-1 || x>=gstr.size()) {
							r.add(null);
							}
						else
							{
							r.add(gstr.get(x));
							}
						}
					t.addList(r);
					}
				t.print(out);
				}
			out.println("<<"+countVariants);
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
			VCFUtils.copyHeaderAndVariantsTo(in, viewer);
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
