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
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceContig;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceGenome;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceGenomeFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import htsjdk.variant.vcf.VCFIterator;



/**

BEGIN_DOC


First described in http://plindenbaum.blogspot.fr/2011/01/my-tool-to-annotate-vcf-files.html.



### Examples




#### Example 1




```
$java -jar  dist/vcfpredictions.jar  \
	-R human_g1k_v37.fasta \
	-k <(curl  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk -F '   ' '{if($2 ~ ".*_.*") next; OFS="        "; gsub(/chr/,"",$2);print;}'  ) \
	I=~/WORK/variations.gatk.annotations.vcf.gz | bgzip > result.vcf.gz

```



#### Annotate 1000 genomes




```
$  curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/vcfpredictions.jar \
  -R  human_g1k_v37.fasta \
  -k <(curl  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk -F '      ' '{if($2 ~ ".*_.*") next; OFS=" "; gsub(/chr/,"",$2);print;}'  )  |\
cut -d '    ' -f 1-8 | grep -A 10 CHROM


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	10583	rs58108140	G	A	100	PASS	AA=.;AC=314;AF=0.14;AFR_AF=0.04;AMR_AF=0.17;AN=2184;ASN_AF=0.13;AVGPOST=0.7707;ERATE=0.0161;EUR_AF=0.21;LDAF=0.2327;PRED=|||||intergenic_variant;RSQ=0.4319;SNPSOURCE=LOWCOV;THETA=0.0046;VT=SNP
1	10611	rs189107123	C	G	100	PASS	AA=.;AC=41;AF=0.02;AFR_AF=0.01;AMR_AF=0.03;AN=2184;ASN_AF=0.01;AVGPOST=0.9330;ERATE=0.0048;EUR_AF=0.02;LDAF=0.0479;PRED=|||||intergenic_variant;RSQ=0.3475;SNPSOURCE=LOWCOV;THETA=0.0077;VT=SNP
1	13302	rs180734498	C	T	100	PASS	AA=.;AC=249;AF=0.11;AFR_AF=0.21;AMR_AF=0.08;AN=2184;ASN_AF=0.02;AVGPOST=0.8895;ERATE=0.0058;EUR_AF=0.14;LDAF=0.1573;PRED=uc010nxq.1|||||intron_variant;RSQ=0.6281;SNPSOURCE=LOWCOV;THETA=0.0048;VT=SNP
1	13327	rs144762171	G	C	100	PASS	AA=.;AC=59;AF=0.03;AFR_AF=0.02;AMR_AF=0.03;AN=2184;ASN_AF=0.02;AVGPOST=0.9698;ERATE=0.0012;EUR_AF=0.04;LDAF=0.0359;PRED=uc010nxq.1|||||intron_variant;RSQ=0.6482;SNPSOURCE=LOWCOV;THETA=0.0204;VT=SNP
1	13957	rs201747181	TC	T	28	PASS	AA=TC;AC=35;AF=0.02;AFR_AF=0.02;AMR_AF=0.02;AN=2184;ASN_AF=0.01;AVGPOST=0.8711;ERATE=0.0065;EUR_AF=0.02;LDAF=0.0788;PRED=uc010nxq.1|||||;RSQ=0.2501;THETA=0.0100;VT=INDEL
1	13980	rs151276478	T	C	100	PASS	AA=.;AC=45;AF=0.02;AFR_AF=0.01;AMR_AF=0.02;AN=2184;ASN_AF=0.02;AVGPOST=0.9221;ERATE=0.0034;EUR_AF=0.02;LDAF=0.0525;PRED=uc010nxq.1|||||3_prime_UTR_variant;RSQ=0.3603;SNPSOURCE=LOWCOV;THETA=0.0139;VT=SNP
1	30923	rs140337953	G	T	100	PASS	AA=T;AC=1584;AF=0.73;AFR_AF=0.48;AMR_AF=0.80;AN=2184;ASN_AF=0.89;AVGPOST=0.7335;ERATE=0.0183;EUR_AF=0.73;LDAF=0.6576;PRED=|||||intergenic_variant;RSQ=0.5481;SNPSOURCE=LOWCOV;THETA=0.0162;VT=SNP
1	46402	rs199681827	C	CTGT	31	PASS	AA=.;AC=8;AF=0.0037;AFR_AF=0.01;AN=2184;ASN_AF=0.0017;AVGPOST=0.8325;ERATE=0.0072;LDAF=0.0903;PRED=|||||intergenic_variant;RSQ=0.0960;THETA=0.0121;VT=INDEL
1	47190	rs200430748	G	GA	192	PASS	AA=G;AC=29;AF=0.01;AFR_AF=0.06;AMR_AF=0.0028;AN=2184;AVGPOST=0.9041;ERATE=0.0041;LDAF=0.0628;PRED=|||||intergenic_variant;RSQ=0.2883;THETA=0.0153;VT=INDEL
1	51476	rs187298206	T	C	100	PASS	AA=C;AC=18;AF=0.01;AFR_AF=0.01;AMR_AF=0.01;AN=2184;ASN_AF=0.01;AVGPOST=0.9819;ERATE=0.0021;EUR_AF=0.01;LDAF=0.0157;PRED=|||||intergenic_variant;RSQ=0.5258;SNPSOURCE=LOWCOV;THETA=0.0103;VT=SNP
```





### See also



 *  http://plindenbaum.blogspot.fr/2013/10/yes-choice-of-transcripts-and-software.html
 *  VcfFilterSequenceOntology
 *  GroupByGene



## History

 *  2013-Dec : moved to a standard arg/argv command line







END_DOC
*/

@Program(name="vcfpredictions",
	description="Basic Variant Effect prediction using ucsc-known gene",
	keywords={"vcf","annotation","prediction"}
	)
public class VCFPredictions extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFPredictions.class).make();
	private enum OutputSyntax {Native,Vep,SnpEff };
	private IntervalTreeMap<List<KnownGene>> knownGenes=null;
	private ReferenceGenome referenceGenome = null;
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-k","--knownGene"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String kgURI =KnownGene.getDefaultUri();

	@Parameter(names={"-soacn","--printsoacn"},description="Print SO:term accession rather than label")
	private boolean print_SO_ACN = false;

	@Parameter(names={"-os","--output-syntax","--syntax"},description="[20180122]output formatting syntax. SnpEff is still not complete.")
	private OutputSyntax outputSyntax = OutputSyntax.Native;

	@Parameter(names={"-R","--reference"},description="[20180122](moved to faidx/DAS). "+ReferenceGenomeFactory.OPT_DESCRIPTION,required=true)
	private String referenceGenomeSource = null;

	/** a sequence with one or more altered base/amino-acid */
	private static class MutedSequence extends DelegateCharSequence
		{
		private final Map<Integer, Character> pos2char=new TreeMap<Integer, Character>();
		MutedSequence(final CharSequence wild)
			{
			super(wild);
			}
		
		void put(final int pos,final char c)
			{
			this.pos2char.put(pos, c);
			}
		
		@Override
		public char charAt(final int i)
			{
			return this.pos2char.getOrDefault(i, getDelegate().charAt(i));
			}		
		}
	
	
	/** a protein, cdna translated with a given genetic code. */
	private static class ProteinCharSequence extends DelegateCharSequence
		{
		private final GeneticCode geneticCode;
		ProteinCharSequence(final GeneticCode geneticCode,final CharSequence cDNA)
			{
			super(cDNA);
			this.geneticCode=geneticCode;
			}
		
		@Override
		public char charAt(int i)
			{
			return geneticCode.translate(
				getDelegate().charAt(i*3+0),
				getDelegate().charAt(i*3+1),
				getDelegate().charAt(i*3+2));
			}	
		
		@Override
		public int length()
			{
			return getDelegate().length()/3;
			}
		}
		
	
	class Annotation
		{
		KnownGene kg;
		Allele alt2;
		String exon_name="";
		String intron_name="";
		Integer position_cds=null;
		String wildAA="";
		String mutAA="";
		String wildCodon="";
		String mutCodon="";
		Integer position_protein=null;
		final Set<SequenceOntologyTree.Term> seqont=new HashSet<>();
		
		
		public String toString()
			{
			final StringBuilder b=new StringBuilder();
			switch(VCFPredictions.this.outputSyntax)
				{
				case Vep:
					{
					//Allele|Feature|Feature_type|Consequence|CDS_position|Protein_position|Amino_acids|Codons
					b.append(alt2==null?"":alt2.getBaseString());
					b.append('|');
					b.append(kg==null?"":kg.getName());
					b.append('|');
					b.append("Transcript");
					b.append('|');
					b.append(this.seqont.stream().map(T->T.getLabel().
							replaceAll("[ ]", "_")).
							collect(Collectors.joining("&")));
					b.append('|');
					b.append(position_cds==null?"":String.valueOf(position_cds+1));
					b.append('|');
					b.append(position_protein==null?"":String.valueOf(position_protein));
					b.append('|');
					b.append(wildAA.isEmpty() && mutAA.isEmpty()?"":wildAA+"/"+mutAA);
					b.append('|');
					b.append(wildCodon.isEmpty() && mutCodon.isEmpty()?"":wildCodon+"/"+mutCodon);
					break;
					}
				case SnpEff:
					{
					b.append(alt2==null?"":alt2.getBaseString());
					b.append('|');
					b.append(this.seqont.stream().map(T->T.getLabel().
							replaceAll("[ ]", "_")).
							collect(Collectors.joining("&")));
					b.append('|');
					b.append("MODIFIER");//impact 
					b.append(kg==null?"":kg.getName());// gene name
					b.append('|');
					b.append(kg==null?"":kg.getName()); //gene id
					b.append('|');//
					b.append("transcript");//Feature_Type
					b.append(kg==null?"":kg.getName());//Feature_ID
					b.append('|');
					b.append("protein_coding");//Transcript_BioType
					b.append('|');
					b.append("");//RAnk
					b.append('|');
					b.append("");//HGVS.c //TODO
					b.append('|');
					b.append("");//HGVS.p 
					b.append('|');
					b.append("");//cDNA.pos / cDNA.length 
					b.append('|');
					b.append("");//CDS.pos / CDS.length 
					b.append('|');
					b.append("");//AA.pos / AA.length 
					b.append('|');
					b.append("");//Distance 
					b.append('|');
					b.append("");//ERRORS / WARNINGS / INFO
					break;
					}
				default:
					{
					b.append(kg==null?"":kg.getName());
					b.append('|');
					b.append(position_cds==null?"":String.valueOf(position_cds));
					b.append('|');
					b.append(position_protein==null?"":String.valueOf(position_protein));
					b.append('|');
					b.append(wildCodon.isEmpty() && mutCodon.isEmpty()?"":wildCodon+"/"+mutCodon);
					b.append('|');
					b.append(wildAA.isEmpty() && mutAA.isEmpty()?"":wildAA+"/"+mutAA);
					b.append("|");
					boolean first=true;
					for(final SequenceOntologyTree.Term t:seqont)
						{
						if(!first) b.append('&');
						first=false;
						if(print_SO_ACN)
							{
							b.append(t.getAcn());
							}
						else
							{
							b.append(t.getLabel().replaceAll("[ ]","_"));
							}
						}
					break;
					}
				}
			return b.toString();
			}
		
		}	
	
	private void loadKnownGenesFromUri() throws IOException
		{
		BufferedReader in=null;
		try {
			if (this.referenceGenome.getDictionary() == null) {
				throw new JvarkitException.FastaDictionaryMissing(this.referenceGenomeSource);
			}
			int n_ignored=0;
			int n_genes = 0;
			this.knownGenes = new IntervalTreeMap<>();
			LOG.info("loading genes");

			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(this.referenceGenome.getDictionary());
			
			in = IOUtils.openURIForBufferedReading(this.kgURI);
			String line;
			final CharSplitter tab = CharSplitter.TAB;
			while ((line = in.readLine()) != null) {
				if (StringUtil.isBlank(line))
					continue;
				final String tokens[] = tab.split(line);
				final KnownGene g = new KnownGene(tokens);
				final String normalizedContig = contigNameConverter.apply(g.getContig());
				if(StringUtil.isBlank(normalizedContig)) {
					++n_ignored;
					continue;
				}
				
				if (this.referenceGenome.getDictionary().getSequence(normalizedContig) == null) {
					++n_ignored;
					continue;
				}
				final int extend_gene_search = 5000; // because we want to set
														// SO:5KB_upstream_variant

				final Interval interval = new Interval(
						normalizedContig,
						Math.max(1, g.getTxStart() + 1 - extend_gene_search),
						g.getTxEnd() + extend_gene_search
						);
				List<KnownGene> L= this.knownGenes.get(interval);
				if(L==null) {
					L=new ArrayList<>(2);
					this.knownGenes.put(interval, L);
				}
				++n_genes;
				L.add(g);
			}
			in.close();
			in = null;
			LOG.info("genes:" + n_genes+" (ignored: "+n_ignored+")");
			}
		finally {
			CloserUtil.close(in);
			}
		}
	private boolean isStop(char c)
		{
		return !Character.isLetter(c);
		}
	
	private  boolean isSimpleBase(final Allele a)
		{
		if(a.isSymbolic() || a.isNoCall()) return false;
		final String s=a.getBaseString();
		if(s.length()!=1) return false;
		switch(s.charAt(0))
			{
			case 'A':case 'a':case 'T':case 't':
			case 'G':case 'g':case 'C':case 'c':
				return true;
			}
		return false;
		}
	
	public static final String TAG="PRED";
	public static enum FORMAT1{TRANSCRIPT,CDSPOS,PROTPOS,CODON,AA,SEQONTOLOGY};
	
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator r, VariantContextWriter w)
		{
		ReferenceContig genomicSequence=null;
		try {
		LOG.info("opening REF:"+this.referenceGenomeSource);
		this.referenceGenome=new ReferenceGenomeFactory().
				open(this.referenceGenomeSource);
		loadKnownGenesFromUri();
		final VCFHeader header= r.getHeader();
		
		final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(this.referenceGenome.getDictionary());
		
		final VCFHeader h2=new VCFHeader(header);
		addMetaData(h2);
		
		switch(this.outputSyntax)
			{
			case Vep:
				{
				h2.addMetaDataLine(new VCFInfoHeaderLine("CSQ",
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Consequence type as predicted by VEP"+
						". Format: Allele|Feature|Feature_type|Consequence|CDS_position|Protein_position|Amino_acids|Codons"
						));
				break;
				}
			case SnpEff:
				{
				h2.addMetaDataLine(new VCFInfoHeaderLine("ANN",
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'"
						));
				break;
				}
			default:
				{
				final StringBuilder format=new StringBuilder();
				for(FORMAT1 f:FORMAT1.values())
					{
					if(format.length()>0) format.append("|"); 
					 format.append(f.name()); 
					}
				
				h2.addMetaDataLine(new VCFInfoHeaderLine(TAG, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,
						"Prediction from "+getClass().getSimpleName()+
						". Format: "+format
						));
				break;
				}
			}
		
        w.writeHeader(h2);

		
		final SequenceOntologyTree soTree=SequenceOntologyTree.getInstance();
		final SequenceOntologyTree.Term so_intron=soTree.getTermByAcn("SO:0001627");
		final SequenceOntologyTree.Term so_exon=soTree.getTermByAcn("SO:0001791");
		final SequenceOntologyTree.Term so_splice_donor=soTree.getTermByAcn("SO:0001575");
		final SequenceOntologyTree.Term so_splice_acceptor=soTree.getTermByAcn("SO:0001574");
		final SequenceOntologyTree.Term so_5_prime_UTR_variant=soTree.getTermByAcn("SO:0001623");
		final SequenceOntologyTree.Term so_3_prime_UTR_variant=soTree.getTermByAcn("SO:0001624");
		final SequenceOntologyTree.Term so_splicing_variant=soTree.getTermByAcn("SO:0001568");
		final SequenceOntologyTree.Term so_stop_lost=soTree.getTermByAcn("SO:0001578");
		final SequenceOntologyTree.Term so_stop_gained=soTree.getTermByAcn("SO:0001587");
		final SequenceOntologyTree.Term so_coding_synonymous=soTree.getTermByAcn("SO:0001819");
		final SequenceOntologyTree.Term so_coding_non_synonymous=soTree.getTermByAcn("SO:0001583");
		final SequenceOntologyTree.Term so_intergenic=soTree.getTermByAcn("SO:0001628");
		final SequenceOntologyTree.Term so_nc_transcript_variant=soTree.getTermByAcn("SO:0001619");
		final SequenceOntologyTree.Term so_non_coding_exon_variant=soTree.getTermByAcn("SO:0001792");
		final SequenceOntologyTree.Term _2KB_upstream_variant=soTree.getTermByAcn("SO:0001636");
		final SequenceOntologyTree.Term _5KB_upstream_variant=soTree.getTermByAcn("SO:0001635");
		final SequenceOntologyTree.Term _5KB_downstream_variant=soTree.getTermByAcn("SO:0001633");
		final SequenceOntologyTree.Term _500bp_downstream_variant=soTree.getTermByAcn("SO:0001634");
		
		
		final ProgressFactory.Watcher<VariantContext> progress=ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
		while(r.hasNext())
			{
			final VariantContext ctx=progress.apply(r.next());
			final String normalizedContig=contigNameConverter.apply(ctx.getContig());
			final List<KnownGene> genes=new ArrayList<>();
			
			
			if(!StringUtil.isBlank(normalizedContig)) {
				for(final List<KnownGene> l2: this.knownGenes.getOverlapping(new Interval(
						normalizedContig,
						ctx.getStart(),
						ctx.getEnd() //1-based
						)))
					{
					genes.addAll(l2);
					}
				}
			final List<Annotation> ctx_annotations=new ArrayList<Annotation>();
			if(genes==null || genes.isEmpty())
				{
				//intergenic
				Annotation a=new Annotation();
				a.seqont.add(so_intergenic);
				ctx_annotations.add(a);
				}
			else
				{
				if(genomicSequence==null || !genomicSequence.hasName(normalizedContig))
					{
					LOG.info("getting genomic Sequence for "+normalizedContig);
					genomicSequence= this.referenceGenome.getContig(normalizedContig);
					if(genomicSequence==null) throw new JvarkitException.ContigNotFoundInDictionary(normalizedContig, this.referenceGenome.getDictionary());
					}
				
				for(final KnownGene gene:genes)
					{
					final GeneticCode geneticCode=GeneticCode.getStandard();
            		
            		
					for(final Allele alt2:ctx.getAlternateAlleles())
						{
						if(alt2.isNoCall()) continue;
						if(alt2.isSymbolic())
							{
							LOG.warn("symbolic allele are not handled... "+alt2.getDisplayString());
							continue;
							}
						if(alt2.isReference()) continue;
						final Annotation annotations=new Annotation();
						annotations.kg=gene;
						annotations.alt2=alt2;
						
						if(gene.isNonCoding())
							{
							annotations.seqont.add(so_nc_transcript_variant);
							continue;
							}
						
						
						ctx_annotations.add(annotations);

		        		StringBuilder wildRNA=null;
		        		ProteinCharSequence wildProt=null;
		        		ProteinCharSequence mutProt=null;
		        		MutedSequence mutRNA=null;
		        		int position_in_cds=-1;
		        		
		        		final int position=ctx.getStart()-1;
		        		if(!String.valueOf(genomicSequence.charAt(position)).equalsIgnoreCase(ctx.getReference().getBaseString()))
		        			{
		        			if(isSimpleBase(ctx.getReference()))
			        			{
			        			LOG.warn("Warning REF!=GENOMIC SEQ!!! at "+position+"/"+ctx.getReference());
			        			}
		        			continue;
		        			}
		        		
		        		if(gene.isPositiveStrand())
		            		{
		        			if(position < gene.getTxStart() - 2000) {
		        				annotations.seqont.add(_5KB_upstream_variant);
		        				}
		        			else if(position < gene.getTxStart()) {
		        				annotations.seqont.add(_2KB_upstream_variant);
		        				}
		        			else if( position >= gene.getTxEnd() + 500) {
		        				annotations.seqont.add(_5KB_downstream_variant);
		        				}
		        			else if( position >= gene.getTxEnd() ) {
		        				annotations.seqont.add(_500bp_downstream_variant);
		        				}
		        			else if(position < gene.getCdsStart())
		            			{
		            			annotations.seqont.add(so_5_prime_UTR_variant);//UTR5
		            			}
		            		else if( gene.getCdsEnd()<= position )
		            			{
		            			annotations.seqont.add(so_3_prime_UTR_variant);
		            			}
		            		else
			            		{
			            		int exon_index=0;
			            		while(exon_index< gene.getExonCount())
			            			{
			            			final KnownGene.Exon exon= gene.getExon(exon_index);
			            			
			            			for(int i= exon.getStart();
			            					i< exon.getEnd();
			            					++i)
			            				{
			            				
			            				if(i==position)
			        						{
			        						annotations.exon_name= exon.getName();
			        						if(exon.isNonCoding())
			            						{
			            						annotations.seqont.add(so_non_coding_exon_variant);
			            						}
			        						}
			            				if(i< gene.getTxStart()) continue;
			            				if(i< gene.getCdsStart()) continue;
			            				if(i>=gene.getCdsEnd()) break;
			        					
			        					if(wildRNA==null)
			        						{
			        						wildRNA=new StringBuilder();
			        						mutRNA=new MutedSequence(wildRNA);
			        						}
			        					
			        					if(i==position)
			        						{
			        						annotations.seqont.add(so_exon);
			        						annotations.exon_name=exon.getName();
			        						position_in_cds=wildRNA.length();
			        						annotations.position_cds= position_in_cds;
			        						//in splicing ?
			        						if(exon.isSplicing(position))
			        							{
			        							
			        							if(exon.isSplicingAcceptor(position))
			        								{
			        								annotations.seqont.add(so_splice_acceptor); //SPLICING_ACCEPTOR
			        								}
			        							else  if(exon.isSplicingDonor(position))
			        								{
			        								annotations.seqont.add(so_splice_donor); // SPLICING_DONOR
			        								}
			        							else //??
			        								{
				        							annotations.seqont.add(so_splicing_variant);
			        								}
			        							}
			        						}
			        					
			            				wildRNA.append(genomicSequence.charAt(i));
			            				
			            				if(i==position && 
			            						isSimpleBase(alt2) && 
			            						isSimpleBase(ctx.getReference()))
			            					{
			            					mutRNA.put(
			            							position_in_cds,
			            							alt2.getBaseString().charAt(0)
			            							);
			            					
			            					}
			            				
			            				if(wildRNA.length()%3==0 && wildRNA.length()>0 && wildProt==null)
				            				{
				            				wildProt=new ProteinCharSequence(geneticCode,wildRNA);
				            				mutProt=new ProteinCharSequence(geneticCode,mutRNA);
				            				}
			            				}
			            			final KnownGene.Intron intron= exon.getNextIntron();
			            			if(intron!=null && intron.contains(position))
			            				{
			            				annotations.intron_name=intron.getName();
			            				annotations.seqont.add(so_intron);
			            				
			            				if(intron.isSplicing(position))
			        						{
			        						if(intron.isSplicingAcceptor(position))
			        							{
			        							annotations.seqont.add(so_splice_acceptor);
			        							}
			        						else if(intron.isSplicingDonor(position))
			        							{
			        							annotations.seqont.add(so_splice_donor);
			        							}
			        						else //???
			        							{
			        							annotations.seqont.add(so_splicing_variant);
			        							}
			        						}
			            				}
			            			++exon_index;
			            			}
			            		}
		            		
		            		
		            		}
		            	else // reverse orientation
		            		{
		            		if(position >= gene.getTxEnd() + 2000) {
		        				annotations.seqont.add(_5KB_upstream_variant);
		        				}
		        			else if(position >= gene.getTxEnd()) {
		        				annotations.seqont.add(_2KB_upstream_variant);
		        				}
		        			else if( position < gene.getTxStart() - 500) {
		        				annotations.seqont.add(_5KB_downstream_variant);
		        				}
		        			else if( position < gene.getTxStart() ) {
		        				annotations.seqont.add(_500bp_downstream_variant);
		        				}
		        			else if(position < gene.getCdsStart())
		            			{
		            			annotations.seqont.add(so_3_prime_UTR_variant);
		            			}
		            		else if( gene.getCdsEnd()<=position )
		            			{
		            			annotations.seqont.add(so_5_prime_UTR_variant);
		            			}
		            		else
			            		{
			            		int exon_index = gene.getExonCount()-1;
			            		while(exon_index >=0)
			            			{

			            			final KnownGene.Exon exon= gene.getExon(exon_index);
			            			
			            			
			            			for(int i= exon.getEnd()-1;
			            				    i>= exon.getStart();
			            				--i)
			            				{
			            				
			            				
			            				if(i==position)
			        						{
			            					annotations.exon_name=exon.getName();
			            					if(exon.isNonCoding())
			            						{
			            						annotations.seqont.add(so_non_coding_exon_variant);
			            						}
			        						}
			            				if(i>= gene.getCdsEnd()) continue;
			            				if(i<  gene.getCdsStart()) break;
			            				
			            				
			            				if(wildRNA==null)
			        						{
			        						wildRNA=new StringBuilder();
			        						mutRNA=new MutedSequence(wildRNA);
			        						}
			            				
			            				if(i==position)
			        						{
			            					annotations.seqont.add(so_exon);
			            					position_in_cds=wildRNA.length();
			        						annotations.position_cds=position_in_cds;
			        						//in splicing ?
			        						if(exon.isSplicing(position))
			        							{		        							
			        							if(exon.isSplicingAcceptor(position))
			        								{
			        								annotations.seqont.add(so_splice_acceptor);
			        								}
			        							else  if(exon.isSplicingDonor(position))
			        								{
			        								annotations.seqont.add(so_splice_donor);
			        								}
			        							else //?
			        								{
			        								annotations.seqont.add(so_splicing_variant);
			        								}
			        							}
			        						
			        						if(isSimpleBase(alt2) &&
			        							isSimpleBase(ctx.getReference()))
				        						{
				        						mutRNA.put(
				        								position_in_cds,
				        								AcidNucleics.complement(alt2.getBaseString().charAt(0))
				        								);
				        						}
			        						}
			            				
			            				wildRNA.append(AcidNucleics.complement(genomicSequence.charAt(i)));
			            				if( wildRNA.length()%3==0 &&
			            					wildRNA.length()>0 &&
			            					wildProt==null)
				            				{
				            				wildProt=new ProteinCharSequence(geneticCode,wildRNA);
				            				mutProt=new ProteinCharSequence(geneticCode,mutRNA);
				            				}
			            				}
			            			final KnownGene.Intron intron= exon.getPrevIntron();
			            			if(intron!=null &&
			            				intron.contains(position))
			            				{
			            				annotations.intron_name=intron.getName();
			            				annotations.seqont.add(so_intron);
			            				
			            				if(intron.isSplicing(position))
			        						{
			        						if(intron.isSplicingAcceptor(position))
			        							{
			        							annotations.seqont.add(so_splice_acceptor);
			        							}
			        						else if(intron.isSplicingDonor(position))
			        							{
			        							annotations.seqont.add(so_splice_donor);
			        							}
			        						else //?	
			        							{
			        							annotations.seqont.add(so_splicing_variant);
			        							}
			        						}
			            				}
			            			--exon_index;
			            			}
			            		}

		            		}//end of if reverse
		        		
		        		
		        		if( isSimpleBase(alt2) &&
		        			isSimpleBase(ctx.getReference()) &&
		        			wildProt!=null &&
		        			mutProt!=null && 
		        			position_in_cds>=0)
			    			{
		            		final int pos_aa=position_in_cds/3;
		            		final int mod= position_in_cds%3;
		            		annotations.wildCodon=(""+
		            			wildRNA.charAt(position_in_cds-mod+0)+
		            			wildRNA.charAt(position_in_cds-mod+1)+
		            			wildRNA.charAt(position_in_cds-mod+2)
		            			);
		            		annotations.mutCodon=(""+
		            			mutRNA.charAt(position_in_cds-mod+0)+
		            			mutRNA.charAt(position_in_cds-mod+1)+
		            			mutRNA.charAt(position_in_cds-mod+2)
		            			);
		            		annotations.position_protein=(pos_aa+1);
		            		annotations.wildAA=String.valueOf(wildProt.charAt(pos_aa));
		            		annotations.mutAA=(String.valueOf(mutProt.charAt(pos_aa)));
		            		
		            		annotations.seqont.remove(so_exon);
		            		
			    			if(isStop(wildProt.charAt(pos_aa)) &&
			    			   !isStop(mutProt.charAt(pos_aa)))
			    				{
			    				annotations.seqont.add(so_stop_lost);
			    				}
			    			else if( !isStop(wildProt.charAt(pos_aa)) &&
			    				 isStop(mutProt.charAt(pos_aa)))
			    				{
			    				annotations.seqont.add(so_stop_gained);
			    				}
			    			else if(wildProt.charAt(pos_aa)==mutProt.charAt(pos_aa))
			    				{
			    				annotations.seqont.add(so_coding_synonymous);
			    				}
			    			else
			    				{
			    				annotations.seqont.add(so_coding_non_synonymous);
			    				}
			    			}
		        		
						}
					}
				}
			
		
			
			final Set<String> info=new HashSet<String>(ctx_annotations.size());
			for(final Annotation a:ctx_annotations)
				{
				info.add(a.toString());
				}
			
			final VariantContextBuilder vb=new VariantContextBuilder(ctx);
			final String thetag;
			switch(this.outputSyntax)
				{
				case Vep : thetag="CSQ"; break;
				case SnpEff : thetag="ANN"; break;
				default: thetag=TAG;break;
				}
			vb.attribute(thetag, info.toArray());
			w.add(vb.make());
			}
		progress.close();
		return RETURN_OK;
		} catch(Exception err ) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(this.referenceGenome);
		}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(StringUtil.isBlank(this.referenceGenomeSource))
			{
			LOG.error("Reference undefined.");
			return -1;
			}
		if(this.kgURI==null || this.kgURI.trim().isEmpty()) 
			{
			LOG.error("knownGene undefined.");
			return -1;
			}
		return doVcfToVcf(args,outputFile);
		}
	
	public static void main(final String[] args)
		{
		new VCFPredictions().instanceMainWithExit(args);
		}
	}
