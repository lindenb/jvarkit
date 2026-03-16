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
package com.github.lindenb.jvarkit.tools.vcf2xls;
/** 
BEGIN_DOC

## Motivation

convert a vcf to something excel. I'm not proud about this. I was asked to do this but please, do not use excel, do not use this tool.

## Example

```
$ java -jar dist/jvarkit.jar vcf2xls src/test/resources/rotavirus_rf.ann.vcf.gz  
#CHROM [1]	POS [1]	END [1]	Allele [1]	Allele [2]	QUAL [1]	FILTER/PASS [1]	INFO/AC [1]	INFO/AN [1]	INFO/ANN [1]	INFO/ANN [2]	INFO/ANN [3]	INFO/BQB [1]	INFO/DP [1]	INFO/DP4 [1]	INFO/DP4 [2]	INFO/DP4 [3]	INFO/DP4 [4]	INFO/HOB [1]	INFO/ICB [1]	INFO/IDV [1]	INFO/IMF [1]	INFO/INDEL [1]	INFO/LOF [1]	INFO/MQ [1]	INFO/MQ0F [1]	INFO/MQB [1]	INFO/MQSB [1]	INFO/RPB [1]	INFO/SGB [1]	INFO/VDB [1]	ANN/Allele [1]	ANN/Allele [2]	ANN/Allele [3]	ANN/Prediction [1]	ANN/Prediction [2]	ANN/Prediction [3]	ANN/GeneId [1]	ANN/GeneId [2]	ANN/GeneId [3]	ANN/GeneName [1]	ANN/GeneName [2]	ANN/GeneName [3]	ANN/FeatureType [1]	ANN/FeatureType [2]	ANN/FeatureType [3]	ANN/FeatureId [1]	ANN/FeatureId [2]	ANN/FeatureId [3]	S1/genotype [1]	S1/type [1]	S1/Allele [1]	S1/Allele [2]	S1/PL [1]	S2/genotype [1]S2/type [1]	S2/Allele [1]	S2/Allele [2]	S2/PL [1]	S3/genotype [1]	S3/type [1]	S3/Allele [1]	S3/Allele [2]	S3/PL [1]	S4/genotype [1]	S4/type [1]	S4/Allele [1]	S4/Allele [2]	S4/PL [1]	S5/genotype [1]S5/type [1]	S5/Allele [1]	S5/Allele [2]	S5/PL [1]
RF01	970	970	A	C	48.67	true	2	10	C|missense_variant|MODERATE|Gene_18_3284|Gene_18_3284|transcript|AAA47319.1|protein_coding|1/1|c.952A>C|p.Lys318Gln|952/3267|952/3267|318/1088||			0.572843	36	19	7	3	5	0.32	0.425					60	0	1	1	0.658863	10.3229	0.693968	C			missense_variant			Gene_18_3284			Gene_18_3284			transcript			AAA47319.1			0/0	HOM_REF	A	A	0,9,47	0/0	HOM_REF	A	A	0,18,73	0/0	HOM_REF	A	A	0,18,730/0	HOM_REF	A	A	0,33,116	1/1	HOM_VAR	C	C	95,24,0
RF02	251	251	A	T	21.29	true	2	10	T|stop_gained|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.235A>T|p.Lys79*|235/2643|235/2643|79/880||	T|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-1371A>T|||||1371|WARNING_TRANSCRIPT_INCOMPLETE	T|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1758A>T|||||1758|WARNING_TRANSCRIPT_NO_START_CODON	1	24	18	0	6	0	0.08	0.235294					60	0	1		0.566154	2.05141	0.0744703	T	T	T	stop_gained	upstream_gene_variant	upstream_gene_variant	UniProtKB/Swiss-Prot:P12472	Gene_1621_1636	UniProtKB/Swiss-Prot:P12472	UniProtKB/Swiss-Prot:P12472	Gene_1621_1636	UniProtKB/Swiss-Prot:P12472	transcript	transcript	transcript	CAA32213.1	CAA32214.1	CAA32215.1	0/0	HOM_REF	A	A	0,15,57	0/1	HET	A	T	31,0,5	0/1	HET	A	T31,0,5	0/0	HOM_REF	A	A	0,9,42	0/0	HOM_REF	A	A	0,24,69
RF02	578	578	G	A	53.0	true	2	10	A|missense_variant|MODERATE|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.562G>A|p.Val188Ile|562/2643|562/2643|188/880||	A|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-1044G>A|||||1044|WARNING_TRANSCRIPT_INCOMPLETE	A|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1431G>A|||||1431|WARNING_TRANSCRIPT_NO_START_CODON	1	57	21	25	8	3	0.32	0.425					60	0	1	1	0.098387	15.148	0.0556216	A	A	A	missense_variant	upstream_gene_variant	upstream_gene_variant	UniProtKB/Swiss-Prot:P12472	Gene_1621_1636	UniProtKB/Swiss-Prot:P12472	UniProtKB/Swiss-Prot:P12472	Gene_1621_1636	UniProtKB/Swiss-Prot:P12472	transcript	transcript	transcript	CAA32213.1	CAA32214.1	CAA32215.1	0/0	HOM_REF	G	G	0,33,122	0/0	HOM_REF	G	G	0,39,135	0/0	HOM_REF	G	G	0,39,135	1/1	HOM_VAR	A	A	100,30,0	0/0	HOM_REF	G	G	0,27,109
RF02	877	877	T	A	3.45	true	1	10	A|missense_variant|MODERATE|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.861T>A|p.Asn287Lys|861/2643|861/2643|287/880||	A|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-745T>A|||||745|WARNING_TRANSCRIPT_INCOMPLETE	A|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1132T>A|||||1132|WARNING_TRANSCRIPT_NO_START_CODON	1	46	19	21	4	2	0.02	0.0439024					60	0	1	1	0.841693	-7.90536	0.479322	A	A	A	missense_variant	upstream_gene_variant	upstream_gene_variant	UniProtKB/Swiss-Prot:P12472	Gene_1621_1636	UniProtKB/Swiss-Prot:P12472	UniProtKB/Swiss-Prot:P12472	Gene_1621_1636	UniProtKB/Swiss-Prot:P12472	transcript	transcript	transcript	CAA32213.1	CAA32214.1	CAA32215.1	0/1	HET	T	A	37,0,50	0/0	HOM_REF	T	T	0,22,116	0/0	HOM_REF	T	T	0,22,116	0/0	HOM_REF	T	T	0,21,94	0/0	HOM_REF	T	T	0,12,62
RF02	1726	1726	T	G	8.23	true	2	10	G|synonymous_variant|LOW|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.1710T>G|p.Thr570Thr|1710/2643|1710/2643|570/880||	G|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-283T>G|||||283|WARNING_TRANSCRIPT_NO_START_CODON	G|downstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.*89T>G|||||89|WARNING_TRANSCRIPT_INCOMPLETE	1	38	18	13	0	7	0.08	0.235294					60	0	1	1	0.870265	-8.29846	0.0568693	G	G	G	synonymous_variant	upstream_gene_variant	downstream_gene_variant	UniProtKB/Swiss-Prot:P12472	UniProtKB/Swiss-Prot:P12472	Gene_1621_1636	UniProtKB/Swiss-Prot:P12472	UniProtKB/Swiss-Prot:P12472	Gene_1621_1636	transcript	transcript	transcript	CAA32213.1	CAA32215.1	CAA32214.1	0/0	HOM_REF	T	T	0,18,83	0/1	HET	T	G	24,0,40	0/1	HET	T	G	24,0,40	0/0	HOM_REF	T	T	0,27,111	0/0	HOM_REF	T	T	0,10,78
```

END_DOC
*/
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
@Program(name="vcf2xls",
description="convert a vcf to something excel. I'm not proud about this. I was asked to do this but please, do not use excel, do not use this tool.",
keywords={"vcf","table","tsv"},
creationDate="20260314",
modificationDate="20260316",
jvarkit_amalgamion =  true,
menu="VCF Manipulation"
)
public class VcfToExcel extends Launcher {
	private static final Logger LOG = Logger.of(VcfToExcel.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-G","--one-genotype-per-line"},description="Print one genotype per line")
	private boolean one_genotype_per_line=false;
	@Parameter(names={"--header"},description="Print header description")
	private boolean print_header =false;
	@Parameter(names={"--empty"},description="Fill Empty cells with that string.")
	private String empty_string="";
	
	private abstract class ColumnHandler {
		VCFHeader vcfHeader;
		AnnPredictionParser annPredictionParser=null;
		
		void initHeader(final VCFHeader vcfHeader) {
			this.vcfHeader = vcfHeader;
			this.annPredictionParser = new AnnPredictionParserFactory(vcfHeader).get();
			}
		
		abstract void beginVariant();
		abstract void endVariant();
		abstract void column(String name,Object v);
		
		private void variantInfo(final VariantContext ctx) {
			column("CHROM",ctx.getContig());
			column("POS",ctx.getStart());
			column("END",ctx.getEnd());
			if(ctx.hasID()) {
				for(String id: CharSplitter.SEMICOLON.split(ctx.getID())) {
					column("ID",id);
					}
				}
			for(Allele a: ctx.getAlleles()) {
				column("Allele",a.getDisplayString());
				}
			if(ctx.hasLog10PError()) {
				column("QUAL",ctx.getPhredScaledQual());
				}
			if(ctx.isNotFiltered()) {
				column("FILTER/PASS",true);
				}
			for(VCFFilterHeaderLine h: this.vcfHeader.getFilterLines() ) {
				if(!ctx.getFilters().contains(h.getID())) continue;
				column("FILTER/"+h.getID(),true);
				}
			for(VCFInfoHeaderLine h: this.vcfHeader.getInfoHeaderLines() ) {
				if(ctx.hasAttribute(h.getID())) {
					for(Object o: ctx.getAttributeAsList(h.getID())) {
						column("INFO/"+h.getID(),o);
						}
					}
				}
			for(final AnnPredictionParser.AnnPrediction pred:this.annPredictionParser.getPredictions(ctx)) {
				column("ANN/Allele",pred.getAllele());
				column("ANN/Prediction",pred.getSOTermsString());
				column("ANN/GeneId", pred.getGeneId());
				column("ANN/GeneName", pred.getGeneName());
				column("ANN/FeatureType", pred.getFeatureType());
				column("ANN/FeatureId", pred.getFeatureId());
				}
			}
		
		private void genotypeInfo(final VariantContext ctx,final String prefix,Genotype g) {
			
			column(prefix+"genotype",g.getAlleles().stream().map(A->A.isNoCall()?".":String.valueOf(ctx.getAlleleIndex(A))).collect(Collectors.joining(g.isPhased()?"|":"/")));
			column(prefix+"type",g.getType().name());
			for(Allele a:g.getAlleles()) {
				column(prefix+"Allele",a.getDisplayString());
				}
			if(g.hasGQ()) {
				column(prefix+"GQ",g.getGQ());
				}
			if(g.hasAD()) {
				column(prefix+"AD",g.getAD());
				}
			if(g.hasDP()) {
				column(prefix+"DP",g.getDP());
				}
			if(g.hasPL()) {
				column(prefix+"PL",g.getPL());
				}
			g.getExtendedAttributes().entrySet().forEach(KV->
				column(prefix+KV.getKey(),KV.getValue())
				);
			}
		
		void visit(final VariantContext ctx) {
			if(ctx.hasGenotypes() && one_genotype_per_line) {
				for(Genotype g: ctx.getGenotypes() ) {
					beginVariant();
					variantInfo(ctx);
					column("sample-name",g.getSampleName());
					genotypeInfo(ctx,"FORMAT/",g);
					endVariant();
					}
				}
			else
				{
				beginVariant();
				variantInfo(ctx);
				for(Genotype g: ctx.getGenotypes() ) {
					genotypeInfo(ctx,g.getSampleName()+"/",g);
					}
				endVariant();
				}
			}
		}
	private  class ColumnCollector extends ColumnHandler {
		private final Counter<String> columnToCount=new Counter<>();
		private Counter<String> current_columnToCount=new Counter<String>();
		protected Set<String> ordered_columns= new LinkedHashSet<>();
		@Override
		void initHeader(VCFHeader vcfHeader) {
			super.initHeader(vcfHeader);
			
			this.ordered_columns.add("CHROM");
			this.ordered_columns.add("POS");
			this.ordered_columns.add("END");
			this.ordered_columns.add("ID");
			this.ordered_columns.add("Allele");
			this.ordered_columns.add("QUAL");
			this.ordered_columns.add("FILTER/PASS");
			for(VCFFilterHeaderLine h: this.vcfHeader.getFilterLines() ) {
				this.ordered_columns.add("FILTER/"+h.getID());
				}
			for(VCFInfoHeaderLine h: this.vcfHeader.getInfoHeaderLines() ) {
				this.ordered_columns.add("INFO/"+h.getID());
				}
			this.ordered_columns.add("ANN/Allele");
			this.ordered_columns.add("ANN/Prediction");
			this.ordered_columns.add("ANN/GeneId");
			this.ordered_columns.add("ANN/GeneName");
			this.ordered_columns.add("ANN/FeatureType");
			this.ordered_columns.add("ANN/FeatureId");
			
			if(one_genotype_per_line) {
				this.ordered_columns.add("sample-name");
				this.ordered_columns.add(("FORMAT/genotype"));
				this.ordered_columns.add(("FORMAT/type"));
				this.ordered_columns.add(("FORMAT/Allele"));

				for(VCFFormatHeaderLine f: this.vcfHeader.getFormatHeaderLines()) {
					this.ordered_columns.add(("FORMAT/"+f.getID()));
					}
				}
			else
				{
				for(String sn: this.vcfHeader.getSampleNamesInOrder()) {
					this.ordered_columns.add((sn+"/genotype"));
					this.ordered_columns.add((sn+"/type"));
					this.ordered_columns.add((sn+"/Allele"));
					for(VCFFormatHeaderLine f: this.vcfHeader.getFormatHeaderLines()) {
						this.ordered_columns.add(sn+"/"+f.getID());
						}
					}
				}
			}
		
		@Override
		void column(final String name,final Object v) {
			this.ordered_columns.add(name);
			this.current_columnToCount.incr(name);
			while(columnToCount.count(name) < current_columnToCount.count(name) ) {
				columnToCount.incr(name);
				}
			}
		@Override
		void beginVariant() {
			current_columnToCount = new Counter<>();
			}
		@Override
		void endVariant() {
			
			}
		}
	
	
	private  class PrintCollector extends ColumnHandler {
		private final ColumnCollector delegate;
		private final PrintWriter pw;
		private final Map<String,List<Object>> col2values = new HashMap<>();
		PrintCollector(final PrintWriter pw,final ColumnCollector delegate) {
			this.pw = pw;
			this.delegate=delegate;
			
			if(print_header) {
				pw.println("## FILTERS " + StringUtils.repeat(20, '#'));
				for(VCFFilterHeaderLine h: delegate.vcfHeader.getFilterLines()) {
					pw.println("##FILTER\t"+h.getID()+"\t"+h.getDescription());
					}
				pw.println("## INFOS : metadata about the variant " + StringUtils.repeat(20, '#'));
				for(VCFInfoHeaderLine h: delegate.vcfHeader.getInfoHeaderLines()) {
					pw.println("##INFO\t"+h.getID()+"\t"+h.getDescription());
					}
				pw.println("## FORMAT : metadata about the genotype " + StringUtils.repeat(20, '#'));
				for(VCFFormatHeaderLine h: delegate.vcfHeader.getFormatHeaderLines()) {
					pw.println("##FORMAT\t"+h.getID()+"\t"+h.getDescription());
					}
				}
			
			boolean first=true;
			for(String label: delegate.ordered_columns) {
				final int expect = (int)delegate.columnToCount.count(label);
				for(int i=0;i< expect;i++) {
					this.pw.print(first?"#":"\t");
					first=false;
					this.pw.print(label+" ["+(i+1)+"]");
					}
				}
			this.pw.println();
			}
		
		void column(String name,Object v) {
			if(!delegate.ordered_columns.contains(name)) throw new IllegalStateException();
			List<Object> L = this.col2values.get(name);
			if(L==null) {
				L=new ArrayList<>();
				this.col2values.put(name, L);
				}
			L.add(v);
			}
		
		void beginVariant() {
			col2values.clear();
			}
		String toString(final Object o) {
			if(o==null) return "N/A";
			if(o.getClass().isArray()) {
				final List<String> L= new ArrayList<String>(Array.getLength(o));
				for(int i=0;i< Array.getLength(o);++i) {
					L.add(toString(Array.get(o, i)));
					}
				return String.join(",", L);
				}
			return String.valueOf(o);
			}
		@Override
		void endVariant() {
			boolean first=true;
			//System.err.println(delegate.ordered_columns);
			//System.err.println(this.col2values.entrySet());
			
			for(String label: this.delegate.ordered_columns) {
				List<Object> L=this.col2values.get(label);
				if(L==null) L=Collections.emptyList();
				for(Object o: L) {
					if(!first) this.pw.print("\t");
					first = false;
					pw.print(toString(o));
					}
				final int expect = (int)delegate.columnToCount.count(label);
				for(int i=L.size();i< expect;i++) {
					if(!first) this.pw.print("\t");
					first = false;
					pw.print(VcfToExcel.this.empty_string);
					}
				}
			this.pw.println();
			}
		}
	
	private void scanVcf(final String path,final ColumnHandler handler) throws IOException {
		try(VCFIterator r= new VCFIteratorBuilder().open(path)) {
			final VCFHeader h=r.getHeader();
			handler.initHeader(h);
			while(r.hasNext()) {
				handler.visit(r.next());
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = oneAndOnlyOneFile(args);
			final ColumnCollector collector = new ColumnCollector();
			scanVcf(input,collector);
			
			try(PrintWriter pw= super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				scanVcf(input,new PrintCollector(pw,collector));
				pw.flush();
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	public static void main(String[] args) {
		new VcfToExcel().instanceMainWithExit(args);
	}
}
