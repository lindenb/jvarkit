/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfcomposite;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.BiPredicate;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.JexlGenotypePredicate;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
/*
BEGIN_DOC

## Input

input is a VCF file annotated with SNPEff or VEP.

## Example

```

$ cat jeter.ped 
X	S1	S2	S3	0	1
X	S2	0	0	0	0
X	S3	0	0	0	0

$ java -jar dist/vcfcomposite.jar -r report.txt -g gene.txt -p jeter.ped src/test/resources/rotavirus_rf.ann.vcf.gz --filter "" | grep -v "##"
[WARN][VepPredictionParser]NO INFO[CSQ] found in header. This VCF was probably NOT annotated with VEP. But it's not a problem if this tool doesn't need to access VEP Annotations.
[INFO][VCFComposite]reading variants and genes
[INFO][VCFComposite]compile per gene
[INFO][VCFComposite]write variants
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF02	877	.	T	A	3.45	PASS	AC=1;AN=10;ANN=A|missense_variant|MODERATE|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.861T>A|p.Asn287Lys|861/2643|861/2643|287/880||,A|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-745T>A|||||745|WARNING_TRANSCRIPT_INCOMPLETE,A|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1132T>A|||||1132|WARNING_TRANSCRIPT_NO_START_CODON;BQB=1;COMPOSITE=gene|Gene_1621_1636|source|ANN_GeneId|pos|1962|ref|TACA|sample|S1,gene|UniProtKB/Swiss-Prot:P12472|source|ANN_GeneId|pos|1962|ref|TACA|sample|S1;DP=46;DP4=19,21,4,2;HOB=0.02;ICB=0.0439024;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.841693;SGB=-7.90536;VDB=0.479322	GT:PL	0/1:37,0,50	0/0:0,22,116	0/0:0,22,116	0/0:0,21,94	0/0:0,12,62
RF02	1962	.	TACA	TA	33.43	PASS	AC=1;AN=10;ANN=TA|frameshift_variant|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.1948_1949delCA|p.His650fs|1948/2643|1948/2643|650/880||,TA|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-45_-44delCA|||||45|WARNING_TRANSCRIPT_NO_START_CODON,TA|downstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.*327_*328delCA|||||327|WARNING_TRANSCRIPT_INCOMPLETE;COMPOSITE=gene|Gene_1621_1636|source|ANN_GeneId|pos|877|ref|T|sample|S1,gene|UniProtKB/Swiss-Prot:P12472|source|ANN_GeneId|pos|877|ref|T|sample|S1;DP=43;DP4=22,11,2,0;HOB=0.02;ICB=0.0439024;IDV=3;IMF=0.3;INDEL;LOF=(UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|2|0.50);MQ=60;MQ0F=0;MQSB=1;SGB=0.810227;VDB=0.373246	GT:PL	0/1:70,0,159	0/0:0,15,225	0/0:0,15,225	0/0:0,27,231	0/0:0,27,168
RF04	887	.	A	G	5.31	PASS	AC=1;AN=10;ANN=G|missense_variant|MODERATE|Gene_9_2339|Gene_9_2339|transcript|AAB07453.1|protein_coding|1/1|c.878A>G|p.Glu293Gly|878/2331|878/2331|293/776||;BQB=1;COMPOSITE=gene|Gene_9_2339|source|ANN_GeneId|pos|1857|ref|CAGA|sample|S1;DP=48;DP4=16,28,3,1;HOB=0.02;ICB=0.0439024;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.90467;SGB=3.91248;VDB=0.811811	GT:PL	0/1:40,0,28	0/0:0,24,98	0/0:0,24,98	0/0:0,33,120	0/0:0,42,134
RF04	1857	.	CAGA	CA	39.47	PASS	AC=1;AN=10;ANN=CA|frameshift_variant|HIGH|Gene_9_2339|Gene_9_2339|transcript|AAB07453.1|protein_coding|1/1|c.1850_1851delGA|p.Arg617fs|1850/2331|1850/2331|617/776||;COMPOSITE=gene|Gene_9_2339|source|ANN_GeneId|pos|887|ref|A|sample|S1;DP=45;DP4=12,21,1,1;HOB=0.02;ICB=0.0439024;IDV=2;IMF=0.166667;INDEL;LOF=(Gene_9_2339|Gene_9_2339|1|1.00);MQ=60;MQ0F=0;MQSB=1;SGB=0.810227;VDB=0.969947	GT:PL	0/1:76,0,152	0/0:0,18,194	0/0:0,18,194	0/0:0,15,166	0/0:0,33,255


$ column -t gene.txt
#CHROM  bed.start  bed.end  gene.key                     gene.label                   gene.source  affected.counts  affected.total  affected.samples
RF02    876        1965     Gene_1621_1636               Gene_1621_1636               ANN_GeneId   1                1               S1
RF04    886        1860     Gene_9_2339                  Gene_9_2339                  ANN_GeneId   1                1               S1
RF02    876        1965     UniProtKB/Swiss-Prot:P12472  UniProtKB/Swiss-Prot:P12472  ANN_GeneId   1                1               S1


$ verticalize report.txt 

>>> 2
$1                       #CHROM : RF02
$2                    bed.start : 876
$3                      bed.end : 1965
$4       count.variants.in.gene : 2
$5                     gene.key : Gene_1621_1636
$6                   gene.label : Gene_1621_1636
$7                  gene.source : ANN_GeneId
$8               variant1.start : 877
$9                 variant1.end : 877
$10                variant1.ref : T
$11                variant1.alt : A
$12               variant1.info : ANN=A|missense_variant|MODERATE|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.861T>A|p.Asn287Lys|861/2643|861/2643|287/880||,A|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-745T>A|||||745|WARNING_TRANSCRIPT_INCOMPLETE,A|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1132T>A|||||1132|WARNING_TRANSCRIPT_NO_START_CODON
$13    variant1.gt[S1].affected : HET
$14  variant1.gt[S2].unaffected : HOM_REF
$15  variant1.gt[S3].unaffected : HOM_REF
$16              variant2.start : 1962
$17                variant2.end : 1965
$18                variant2.ref : TACA
$19                variant2.alt : TA
$20               variant2.info : ANN=TA|frameshift_variant|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.1948_1949delCA|p.His650fs|1948/2643|1948/2643|650/880||,TA|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-45_-44delCA|||||45|WARNING_TRANSCRIPT_NO_START_CODON,TA|downstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.*327_*328delCA|||||327|WARNING_TRANSCRIPT_INCOMPLETE
$21    variant2.gt[S1].affected : HET
$22  variant2.gt[S2].unaffected : HOM_REF
$23  variant2.gt[S3].unaffected : HOM_REF
<<< 2

>>> 3
$1                       #CHROM : RF04
$2                    bed.start : 886
$3                      bed.end : 1860
$4       count.variants.in.gene : 2
$5                     gene.key : Gene_9_2339
$6                   gene.label : Gene_9_2339
$7                  gene.source : ANN_GeneId
$8               variant1.start : 887
$9                 variant1.end : 887
$10                variant1.ref : A
$11                variant1.alt : G
$12               variant1.info : ANN=G|missense_variant|MODERATE|Gene_9_2339|Gene_9_2339|transcript|AAB07453.1|protein_coding|1/1|c.878A>G|p.Glu293Gly|878/2331|878/2331|293/776||
$13    variant1.gt[S1].affected : HET
$14  variant1.gt[S2].unaffected : HOM_REF
$15  variant1.gt[S3].unaffected : HOM_REF
$16              variant2.start : 1857
$17                variant2.end : 1860
$18                variant2.ref : CAGA
$19                variant2.alt : CA
$20               variant2.info : ANN=CA|frameshift_variant|HIGH|Gene_9_2339|Gene_9_2339|transcript|AAB07453.1|protein_coding|1/1|c.1850_1851delGA|p.Arg617fs|1850/2331|1850/2331|617/776||
$21    variant2.gt[S1].affected : HET
$22  variant2.gt[S2].unaffected : HOM_REF
$23  variant2.gt[S3].unaffected : HOM_REF
<<< 3

>>> 4
$1                       #CHROM : RF02
$2                    bed.start : 876
$3                      bed.end : 1965
$4       count.variants.in.gene : 2
$5                     gene.key : UniProtKB/Swiss-Prot:P12472
$6                   gene.label : UniProtKB/Swiss-Prot:P12472
$7                  gene.source : ANN_GeneId
$8               variant1.start : 877
$9                 variant1.end : 877
$10                variant1.ref : T
$11                variant1.alt : A
$12               variant1.info : ANN=A|missense_variant|MODERATE|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.861T>A|p.Asn287Lys|861/2643|861/2643|287/880||,A|upstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.-745T>A|||||745|WARNING_TRANSCRIPT_INCOMPLETE,A|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-1132T>A|||||1132|WARNING_TRANSCRIPT_NO_START_CODON
$13    variant1.gt[S1].affected : HET
$14  variant1.gt[S2].unaffected : HOM_REF
$15  variant1.gt[S3].unaffected : HOM_REF
$16              variant2.start : 1962
$17                variant2.end : 1965
$18                variant2.ref : TACA
$19                variant2.alt : TA
$20               variant2.info : ANN=TA|frameshift_variant|HIGH|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32213.1|protein_coding|1/1|c.1948_1949delCA|p.His650fs|1948/2643|1948/2643|650/880||,TA|upstream_gene_variant|MODIFIER|UniProtKB/Swiss-Prot:P12472|UniProtKB/Swiss-Prot:P12472|transcript|CAA32215.1|protein_coding||c.-45_-44delCA|||||45|WARNING_TRANSCRIPT_NO_START_CODON,TA|downstream_gene_variant|MODIFIER|Gene_1621_1636|Gene_1621_1636|transcript|CAA32214.1|protein_coding||c.*327_*328delCA|||||327|WARNING_TRANSCRIPT_INCOMPLETE
$21    variant2.gt[S1].affected : HET
$22  variant2.gt[S2].unaffected : HOM_REF
$23  variant2.gt[S3].unaffected : HOM_REF
<<< 4

```



END_DOC
*/
@Program(name="vcfcomposite",
	description="(in developpement) Finds Variants involved in a Het Compound Disease",
	keywords={"vcf","disease","annotation","pedigree","haplotype"},
	creationDate= "20170331",
	modificationDate = "20200210"
)
public class VCFComposite extends Launcher {
	private static final String INFO_TAG="COMPOSITE";
	private static final Logger LOG= Logger.build(VCFComposite.class).make();
	
	private enum Selection {any,all};
	
	@Parameter(names={"-p","-ped","--pedigree"},description=PedigreeParser.OPT_DESC,required=true)
	private Path pedigreeFile=null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-vf","--variant-filter"},description=JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantJexl = JexlVariantPredicate.create("");
	@Parameter(names={"-gf","--genotype-filter"},description=JexlGenotypePredicate.PARAMETER_DESCRIPTION,converter=JexlGenotypePredicate.Converter.class)
	private BiPredicate<VariantContext,Genotype> genotypeFilter = JexlGenotypePredicate.create("");
	@Parameter(names={"-max","--max","--max-variants"},description="[20180718] Max variants per gene. If different than -1, used to set a maximum number of variants per gene; The idea is to filter out the gene having a large number of variants.")
	private int max_number_of_variant_per_gene = -1;
	@Parameter(names={"--filter"},description="[20180718] set FILTER for the variants that are not part of a composite mutation. Blank = filter out non-composites")
	private String filterNonCompositeTag = "NOT_COMPOSITE";
	@Parameter(names={"-l","--list"},description= "list all available gene extractors", help=true)
	private boolean list_extractors = false;
	@Parameter(names={"-e","-E","--extractors"},description=GeneExtractorFactory.OPT_DESC)
	private String extractorsNames="ANN/GeneId VEP/GeneId";
	@Parameter(names={"-r","--report"},description="Optional tabular text report for pairs of variants")
	private Path reportPath = null;
	@Parameter(names={"-g","--genes"},description="Optional tabular text report for genes")
	private Path geneReportPath = null;
	@Parameter(names={"-s","--select-variant"},description="How to select affected samples for *one* variant. This variant will be later challenged with another variant of the same gene. "
			+ "any: at least one affected sample must be HET for the variant "
			+ "all: all affected samples must be HET for the variant.")
	private Selection affected_selection = Selection.any;
	@Parameter(names={"-u","--one-unaffected-het"},description="There must be at least one **Unaffected** sample with a HET Genotype for a given variant. ")
	private boolean one_unaffected_het = false;
	@Parameter(names={"-a","--select-pair"},description="How to select a pair of variants : "
			+ "any: at least one affected sample must be carried by variant1 and variant2 "
			+ "all: all affected samples be carried by variant1 and variant2")
	private Selection pair_selection = Selection.any;

	
	//@Parameter(names={"-m","--model"},description="Model type",required=true)
	//private Type modelType=null;
	//@Parameter(names={"-models"},description="List the available models and exits",help=true)
	//private boolean listModels=false;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();

	
	private AbstractVCFCodec vcfDecoder = null;
	private VCFEncoder vcfEncoder = null;
	private List<Sample> affectedSamples = new ArrayList<>();
	private List<Sample> unaffectedSamples = new ArrayList<>();
	private PrintWriter reportWriter =  null;
	private PrintWriter reportGeneWriter =  null;
	
	private class VariantLine
		{
		final long id;
		VariantContext ctx;
		VariantLine(long id,VariantContext ctx) {
			this.id = id;
			this.ctx = ctx;
			}
		}
	
	private class VariantLineCodec extends AbstractDataCodec<VariantLine>
		{
		@Override
		public VariantLine decode(final DataInputStream dis) throws IOException {
			try {
				long n = dis.readLong();
				final VariantContext ctx = VCFComposite.this.vcfDecoder.decode(IOUtils.readString(dis));
				return new VariantLine(n,ctx);
				} 
			catch(final EOFException err)
				{
				return null;
				}
			}

		@Override
		public void encode(DataOutputStream dos, VariantLine object) throws IOException {
			dos.writeLong(object.id);
			IOUtils.writeString(dos,VCFComposite.this.vcfEncoder.encode(object.ctx));
			}

		@Override
		public VariantLineCodec clone() {
			return new VariantLineCodec();
			}

		}
	
	private class GeneIdentifier
		implements Comparable<GeneIdentifier>
		{
		String geneName;
		String label;
		String contig="";//not null please
		String source;
		GeneIdentifier() {
			}
		GeneIdentifier(final String geneName,final String label,final String source) {
			this.geneName = geneName;
			this.label = StringUtil.isBlank(label)?".":label;
			this.source = source;
			}

		@Override
		public int hashCode() {
			return this.geneName.hashCode()*31 + this.source.hashCode();
			}
		@Override
		public boolean equals(final Object obj) {
			return this.compareTo(GeneIdentifier.class.cast(obj))==0;
			}
		@Override
		public int compareTo(final GeneIdentifier o) {
			int i = this.geneName.compareTo(o.geneName);
			if(i!=0) return i;
			i = this.contig.compareTo(o.contig);
			if(i!=0) return i;
			i = this.source.compareTo(o.source);
			return i;
			}
		@Override
		public String toString() {
			return geneName+"["+this.label+"]("+source+") on "+contig;
			}
		}

	private class GeneIdentifierCodec
		extends AbstractDataCodec<GeneIdentifier>
		{
		@Override
		public GeneIdentifier decode(final DataInputStream in) throws IOException {
			String gn;
			try {
				gn = in.readUTF();
				} 
			catch(final EOFException err) 
				{
				return null;
				}
			final GeneIdentifier gi = new GeneIdentifier();
			gi.geneName = gn;
			gi.label = in.readUTF();
			gi.contig = in.readUTF();
			gi.source = in.readUTF();
			return gi;
			}
		
		@Override
		public void encode( final DataOutputStream dos, final GeneIdentifier object) throws IOException {
			dos.writeUTF(object.geneName);
			dos.writeUTF(object.label);
			dos.writeUTF(object.contig);
			dos.writeUTF(object.source);
			}
		
		@Override
		public GeneIdentifierCodec clone() {
			return new GeneIdentifierCodec();
			}
		}
	
	private class GeneAndVariant
		{
		final GeneIdentifier gene;
		final VariantLine variant;
				
		public GeneAndVariant(final GeneIdentifier geneIdentifier, final VariantLine variant) {
			this.gene = geneIdentifier;
			if(this.gene==null) throw new IllegalArgumentException("gene is null");
			this.variant = variant;
			if(this.variant==null) throw new IllegalArgumentException("variant is null");
			}
		
		int compareGene(final GeneAndVariant other) {
			return this.gene.compareTo(other.gene);
			}
		
		int compareGeneThenIndex(final GeneAndVariant other) {
			final int i = compareGene(other);
			if(i!=0) return i;
			return Long.compare(this.variant.id, other.variant.id);
			}
		}
	
	private class GeneAndVariantCodec
		extends AbstractDataCodec<GeneAndVariant>
		{
		final GeneIdentifierCodec geneCodec = new GeneIdentifierCodec();
		final VariantLineCodec vcCodec = new VariantLineCodec();
		@Override
		public GeneAndVariant decode(final DataInputStream dis) throws IOException {
			final GeneIdentifier gi;

			try {
				gi  = geneCodec.decode(dis);
				if(gi==null) return null;
				}
			catch(final EOFException err)
				{
				return null;
				}
			final VariantLine vl = vcCodec.decode(dis);
			return new GeneAndVariant(gi, vl);

			}
		@Override
		public void encode(DataOutputStream dos, GeneAndVariant object) throws IOException {
			geneCodec.encode(dos,object.gene);
			vcCodec.encode(dos,object.variant);
			}
		@Override
		public GeneAndVariantCodec clone() {
			return new GeneAndVariantCodec();
			}
		}
	
	private boolean isGenotypeForAffected(final Genotype gt) {
		return gt!=null && gt.isHet();
	}
		
	private final boolean acceptVariant(final VariantContext ctx) {
		/* discard if a control is HOM_VAR */
		if( this.unaffectedSamples.stream().
			map(S->ctx.getGenotype(S.getId())).
			filter(G->G!=null).
			anyMatch(G->G.isHomVar())) return false;
		
		/* at least one het in the unaffected */
		if(this.one_unaffected_het) {
			if( this.unaffectedSamples.stream().
					map(S->ctx.getGenotype(S.getId())).
					filter(G->G!=null).
					noneMatch(G->G.isHet())) return false;
			}
		
		/* case must be het */
		switch(this.affected_selection)
			{
			case any:
				if( this.affectedSamples.stream().
						map(S->ctx.getGenotype(S.getId())).
						anyMatch(G->isGenotypeForAffected(G))) return true;
				break;
			case all:
				if( this.affectedSamples.stream().
						map(S->ctx.getGenotype(S.getId())).
						allMatch(G->isGenotypeForAffected(G))) return true;
				break;
			}
		
		
		return false;
		}
	
		
	
		
	private void scan(
				final GeneIdentifier geneKey,
				final List<VariantLine> variants)
		{
		int pair_index = 0;
		if(variants.size()<2) {
			return;
			}
		// test if this gene is a big pool of variant, false positive. e.g: highly variables genes.
		if(max_number_of_variant_per_gene>=0) {
			if(variants.size() > max_number_of_variant_per_gene) {
				LOG.warn("Too many variants "+variants.size()+" for "+geneKey);
				return;
			}
		}
		
		final Set<Sample> samples_in_gene = new HashSet<>(this.affectedSamples.size());
		
		// search for the one snp
		for(int x=0;x+1< variants.size();++x)
			{
			final VariantContext vcx = variants.get(x).ctx;
			for(int y=x+1;y< variants.size();++y)
				{
				final VariantContext vcy = variants.get(y).ctx;
				
				final Set<Sample> matching_affected_samples = new HashSet<>(this.affectedSamples.size());
				// loop over affected samples
				for(final Sample affected: this.affectedSamples) {
					final Genotype gcx = vcx.getGenotype(affected.getId());
					// child variant  n. y  must be HOM_VAR or HET
					if(gcx==null || !isGenotypeForAffected(gcx)) continue;
					// filtered ?
					if(!this.genotypeFilter.test(vcx,gcx)) continue;

					final Genotype gcy = vcy.getGenotype(affected.getId());
					// child variant n. y must be HOM_VAR or HET
					if(gcy==null || !isGenotypeForAffected(gcy)) continue;
					// filtered ?
					if(!this.genotypeFilter.test(vcy,gcy)) continue;

					
					boolean unaffected_are_ok=true;
					//check unaffected indididual don't have same haplotype
					for(final Sample unaffected: this.unaffectedSamples) {
						final Genotype gux = vcx.getGenotype(unaffected.getId());
						if(gux!=null && gux.isHomVar()) {
							unaffected_are_ok=false;
							break;
							}
						
						final Genotype guy = vcy.getGenotype(unaffected.getId());
						
						if(guy!=null && guy.isHomVar()) {
							unaffected_are_ok=false;
							break;
							}
						
						if(gux!=null && guy!=null &&
							gux.sameGenotype(gcx, true) &&
							guy.sameGenotype(gcy, true)
							)
							{
							unaffected_are_ok=false;
							break;
							}
						}
					if(!unaffected_are_ok) continue;
					
					matching_affected_samples.add(affected);
					}
				
				boolean affected_are_ok;
				
				switch(this.pair_selection)
					{
					case all: affected_are_ok = matching_affected_samples.size()==this.affectedSamples.size() ; break;
					case any: affected_are_ok = !matching_affected_samples.isEmpty();break;
					default: throw new IllegalStateException();
					}
					
				if(affected_are_ok)
					{
					samples_in_gene.addAll(matching_affected_samples);
					
					if(this.reportPath!=null) {
						this.reportWriter.print(geneKey.contig);
						this.reportWriter.print('\t');
						this.reportWriter.print(Math.min(vcx.getStart(),vcy.getStart())-1);
						this.reportWriter.print('\t');
						this.reportWriter.print(Math.max(vcx.getEnd(),vcy.getEnd()));
						this.reportWriter.print('\t');
						this.reportWriter.print(++pair_index);
						this.reportWriter.print('\t');
						this.reportWriter.print(geneKey.geneName);
						this.reportWriter.print('\t');
						this.reportWriter.print(geneKey.label);
						this.reportWriter.print('\t');
						this.reportWriter.print(geneKey.source);
						
						for(int side=0;side<2;++side) 
							{
							final VariantContext vc=(side==0?vcx:vcy);
							this.reportWriter.print('\t');
							this.reportWriter.print(vc.getStart());
							this.reportWriter.print('\t');
							this.reportWriter.print(vc.getEnd());
							this.reportWriter.print('\t');
							this.reportWriter.print(vc.getReference().getDisplayString());
							this.reportWriter.print('\t');
							this.reportWriter.print(vc.getAlternateAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(",")));
							this.reportWriter.print('\t');
							this.reportWriter.print(vc.getAttributes().keySet().
									stream().
									filter(K->K.equals(AnnPredictionParser.getDefaultTag()) || K.equals(VepPredictionParser.getDefaultTag())).
									map(K->K+"="+String.join(",",vc.getAttributeAsStringList(K, ""))).
									collect(Collectors.joining(";"))
									);
							for(final Sample sn:this.affectedSamples) {
								this.reportWriter.print('\t');
								this.reportWriter.print(vc.getGenotype(sn.getId()).getType().name());
								}
							for(final Sample sn:this.unaffectedSamples) {
								this.reportWriter.print('\t');
								this.reportWriter.print(vc.getGenotype(sn.getId()).getType().name());
								}
							}
						this.reportWriter.println();
						}
					
					for(int side=0;side<2;++side) 
						{
						final VariantContext vc=(side==0?vcx:vcy);
						final Set<String> set = getAnnotationsForVariant(vc);
						final VariantContextBuilder vcb = new VariantContextBuilder(vc);
						
						
						final StringBuilder sb=new StringBuilder();
						sb.append("gene|").append(geneKey.geneName);
						sb.append("|name|").append(geneKey.label);
						sb.append("|source|").append(geneKey.source);
						if(side==1) {
							sb.append("|pos|").append(vcx.getStart());
							sb.append("|ref|").append(vcx.getReference().getDisplayString());
							}
						else {
							sb.append("|pos|").append(vcy.getStart());
							sb.append("|ref|").append(vcy.getReference().getDisplayString());
							}
						sb.append("|samples|").append(matching_affected_samples.stream().map(S->S.getId()).collect(Collectors.joining("&")));
						set.add(sb.toString());
						set.remove("");
						vcb.attribute(INFO_TAG, new ArrayList<>(set));
						
						if(side==0)
							{
							variants.get(x).ctx = vcb.make();
							}
						else
							{
							variants.get(y).ctx = vcb.make();
							}
						}
					}
				}
			}
		if(!samples_in_gene.isEmpty() && this.geneReportPath!=null)
			{
			this.reportGeneWriter.print(geneKey.contig);
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(variants.stream().mapToInt(V->V.ctx.getStart()-1).min().orElse(-1));
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(variants.stream().mapToInt(V->V.ctx.getEnd()).max().orElse(-1));
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(geneKey.geneName);
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(geneKey.label);
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(geneKey.source);
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(variants.size());
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(samples_in_gene.size());
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(this.affectedSamples.size());
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print(samples_in_gene.stream().map(S->S.getId()).collect(Collectors.joining(",")));
			this.reportGeneWriter.println();
			}
		}
		

	
	private Set<String> getAnnotationsForVariant(VariantContext vc) {
		return new HashSet<>(vc.getAttributeAsStringList(INFO_TAG, ""));
		}
	
	protected int doVcfToVcf(final String inputName,
		final VCFIterator iterin,
		final VariantContextWriter out) {

		final VCFHeader header = iterin.getHeader();
		if(!header.hasGenotypingData()) {
			LOG.error("No genotypes in "+inputName);
			return -1;
			}
		final GeneExtractorFactory geneExtractorFactory = new GeneExtractorFactory(header);
		final List<GeneExtractorFactory.GeneExtractor> extractors = geneExtractorFactory.parse(this.extractorsNames);
		if(extractors.isEmpty())  {
			LOG.error("no gene extractor found/defined.");
			return -1;
			}	
		
		
		
		final Pedigree pedigree;
		try {
			final Set<String> sampleNames = new HashSet<>(header.getSampleNamesInOrder());
			final PedigreeParser pedParser = new PedigreeParser();
			pedigree = pedParser.parse(this.pedigreeFile);
			if(pedigree==null || pedigree.isEmpty()) {
				LOG.error("pedigree missing/empty");
				return -1;
				}
			this.affectedSamples.addAll(pedigree.getAffectedSamples());
			this.affectedSamples.removeIf(S->!sampleNames.contains(S.getId()));
			if(this.affectedSamples.isEmpty()) {
				LOG.error("No Affected sample in pedigree. "  + this.pedigreeFile+"/"+inputName);
				return -1;
				}
			this.unaffectedSamples.addAll(pedigree.getUnaffectedSamples());
			this.unaffectedSamples.removeIf(S->!sampleNames.contains(S.getId()));
			if(pedigree.getUnaffectedSamples().isEmpty()) {
				LOG.error("No Unaffected sample in " + this.pedigreeFile+"/"+inputName);
				return -1;
				}
			
			}
		catch(final IOException err)
			{
			throw new RuntimeIOException(err);
			}
			
		
		
		
		
		
		header.addMetaDataLine(new VCFInfoHeaderLine(INFO_TAG, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"Variant of VCFComposite"));
		if(!StringUtils.isBlank(this.filterNonCompositeTag)) {
			header.addMetaDataLine(new VCFFilterHeaderLine(this.filterNonCompositeTag, "Not a Variant fir VCFComposite"));
			}
		final SAMSequenceDictionary dict = header.getSequenceDictionary();
		final Comparator<String> contigCmp;
		if(dict==null || dict.isEmpty())
			{
			contigCmp = (A,B)->A.compareTo(B); 
			}
		else
			{
			contigCmp = new ContigDictComparator(dict);
			}
		
		final Comparator<VariantContext> ctxComparator = (V1,V2)->{
			int i = contigCmp.compare(V1.getContig(), V2.getContig());
			if(i!=0) return i;
			i = Integer.compare(V1.getStart(), V2.getStart());
			if(i!=0) return i;
			return V1.getReference().compareTo(V2.getReference());
			};
		
		final Comparator<VariantLine> variantLineComparator = (V1,V2)->{
			final int i = ctxComparator.compare(V1.ctx,V2.ctx);
			if(i!=0) return i;
			return Long.compare(V1.id, V2.id);
			};


		long ID_GENERATOR = 0L;
		this.vcfDecoder = VCFUtils.createDefaultVCFCodec();
		this.vcfDecoder.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
		this.vcfEncoder = new VCFEncoder(header, false, true);
		SortingCollection<GeneAndVariant> sorting=null;
		SortingCollection<VariantLine> outputSorter = null;

		try
			{
			LOG.info("reading variants and genes");
			/* Gene and variant sorter */
			sorting = SortingCollection.newInstance(GeneAndVariant.class,
					new GeneAndVariantCodec(),
					GeneAndVariant::compareGeneThenIndex,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorting.setDestructiveIteration(true);
			/* Variant sorter */
			outputSorter = SortingCollection.newInstance(
					VariantLine.class,
					new VariantLineCodec(),
					variantLineComparator,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			outputSorter.setDestructiveIteration(true);
			
			/* read input */
			while(iterin.hasNext()) {
				final VariantContext ctx = iterin.next();
				final VariantLine variantLine = new VariantLine(++ID_GENERATOR,ctx);
				if(!this.variantJexl.test(ctx))
					{
					outputSorter.add(variantLine);
					continue;
					}
				
				if(!acceptVariant(ctx))
					{
					outputSorter.add(variantLine);
					continue;
					}
				
				final Set<GeneIdentifier> geneKeys = new HashSet<>();
				
				extractors.
						stream().
						map(EX->EX.apply(ctx)).
						flatMap(H->H.keySet().stream()).forEach(KG->
						{
						geneKeys.add(new GeneIdentifier(KG.getKey(),KG.getGene(),KG.getMethod().replace('/', '_')));
						});
				
				if(geneKeys.isEmpty()) {
					outputSorter.add(variantLine);
					continue;
					}
				
				for(final GeneIdentifier gk:geneKeys) 
					{
					final GeneAndVariant gav = new GeneAndVariant(gk,variantLine);
					gav.gene.contig = ctx.getContig();
					sorting.add(gav);
					}
				}
			sorting.doneAdding();
			
			LOG.info("compile per gene");
			
			
			this.reportGeneWriter = (this.geneReportPath==null?new PrintWriter(new NullOuputStream()):IOUtils.openPathForPrintWriter(this.geneReportPath));
			this.reportGeneWriter.print("#CHROM");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("bed.start");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("bed.end");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("gene.key");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("gene.label");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("gene.source");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("count.variants");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("affected.counts");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("affected.total");
			this.reportGeneWriter.print('\t');
			this.reportGeneWriter.print("affected.samples");
			this.reportGeneWriter.println();
			
			this.reportWriter = (this.reportPath==null?new PrintWriter(new NullOuputStream()):IOUtils.openPathForPrintWriter(this.reportPath));

			this.reportWriter.print("#CHROM");
			this.reportWriter.print('\t');
			this.reportWriter.print("bed.start");
			this.reportWriter.print('\t');
			this.reportWriter.print("bed.end");
			this.reportWriter.print('\t');
			this.reportWriter.print("gene.index");
			this.reportWriter.print('\t');
			this.reportWriter.print("gene.key");
			this.reportWriter.print('\t');
			this.reportWriter.print("gene.label");
			this.reportWriter.print('\t');
			this.reportWriter.print("gene.source");
			for(int side=0;side<2;++side) 
				{
				this.reportWriter.print('\t');
				final String prefix="variant"+(side+1)+".";
				this.reportWriter.print(prefix+"start");
				this.reportWriter.print('\t');
				this.reportWriter.print(prefix+"end");
				this.reportWriter.print('\t');
				this.reportWriter.print(prefix+"ref");
				this.reportWriter.print('\t');
				this.reportWriter.print(prefix+"alt");
				this.reportWriter.print('\t');
				this.reportWriter.print(prefix+"info");
				for(final Sample sn:this.affectedSamples) {
					this.reportWriter.print('\t');
					this.reportWriter.print(prefix+"gt["+sn.getId()+"].affected");
					}
				for(final Sample sn:this.unaffectedSamples) {
					this.reportWriter.print('\t');
					this.reportWriter.print(prefix+"gt["+sn.getId()+"].unaffected");
					}
				}
			this.reportWriter.println();
			
			
			//compile data
			CloseableIterator<GeneAndVariant> iter2=sorting.iterator();
			EqualRangeIterator<GeneAndVariant> eqiter= new EqualRangeIterator<>(iter2,(A,B)->A.gene.compareTo(B.gene));
			while(eqiter.hasNext())
				{
				final List<GeneAndVariant> variants=eqiter.next();
				scan(
					variants.get(0).gene,
					variants.stream().
						map(L->L.variant).
						collect(Collectors.toList())
					);
				for(final GeneAndVariant ga:variants) outputSorter.add(ga.variant);
				}
			eqiter.close();
			iter2.close();
			sorting.cleanup();
			//
			this.reportWriter.flush();
			this.reportWriter.close();
			this.reportGeneWriter.flush();
			this.reportGeneWriter.close();
			
			LOG.info("write variants");
			CloseableIterator<VariantLine> iter1 = outputSorter.iterator();
			EqualRangeIterator<VariantLine > eqiter1 = new EqualRangeIterator<>(iter1,variantLineComparator);
			out.writeHeader(header);
			while(eqiter1.hasNext())
				{
				final List<VariantLine> array = eqiter1.next();
				final VariantContext firstCtx = array.get(0).ctx;
				final Set<String> set= getAnnotationsForVariant(firstCtx);
				final VariantContext outCtx;
				final VariantContextBuilder vcb = new VariantContextBuilder(firstCtx);
				for(int y=1;y<array.size();++y) {
					set.addAll(getAnnotationsForVariant(array.get(y).ctx));
					}
				if(set.isEmpty())
					{
					if(StringUtils.isBlank(this.filterNonCompositeTag))
						{
						//ignore
						continue;
						}
					else
						{
						vcb.filter(this.filterNonCompositeTag);
						}
					}
				else
					{
					if(!firstCtx.isFiltered())
						{
						vcb.passFilters();
						}
					vcb.attribute(INFO_TAG, new ArrayList<>(set));
					}

				outCtx = vcb.make();
					
				out.add(outCtx);
				}
			outputSorter.cleanup();
			eqiter1.close();
			iter1.close();
			
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		
		if(this.list_extractors) {
			for(final String en: GeneExtractorFactory.getExtractorNames()) {
				System.out.println(en);
				}
			return 0;
			}
		
		try {
			return doVcfToVcfPath(args,this.writingVariantsDelegate,this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		}
	public static void main(final String[] args) {
		new VCFComposite().instanceMainWithExit(args);
		}
	
	}
