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
package com.github.lindenb.jvarkit.tools.vcfannot;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.Set;
import java.util.WeakHashMap;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.ExonOrIntron;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.PeptideSequence;
import com.github.lindenb.jvarkit.util.bio.structure.RNASequence;
import com.github.lindenb.jvarkit.util.bio.structure.RNASequenceFactory;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.bio.structure.UTR;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;



/**

BEGIN_DOC

First described in http://plindenbaum.blogspot.fr/2011/01/my-tool-to-annotate-vcf-files.html.

### Examples


#### Example 1


```
$java -jar  dist/vcfpredictions.jar  \
	-R human_g1k_v37.fasta \
	--gtf Homo_sapiens.GRCh37.75.gtf \
	input.vcf.gz | bgzip > result.vcf.gz

```

#### Annotate 1000 genomes


```
$  curl "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz" |\
gunzip -c |\
java -jar dist/vcfpredictions.jar \
  -R  human_g1k_v37.fasta \
  --gtf Homo_sapiens.GRCh37.75.gtf  |\
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


END_DOC
*/

@Program(name="vcfpredictions",
	description="Basic Variant Effect prediction using gtf",
	keywords={"vcf","annotation","prediction","gtf"},
	creationDate="20160318",
	modificationDate="20200701"
	)
public class VCFPredictions extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VCFPredictions.class).make();
	private enum OutputSyntax {Native,Vep,SnpEff };
	private IntervalTreeMap<List<Transcript>> transcriptTreeMap = null;
	private ReferenceSequenceFile referenceGenome = null;
	private final WeakHashMap<String, RNASequence> transcriptId2cdna = new WeakHashMap<>();
	private GenomicSequence genomicSequence = null;

	@Parameter(names={"-k","-g","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfPath = null;

	@Parameter(names={"-os","--output-syntax","--syntax"},description="Output formatting syntax.")
	private OutputSyntax outputSyntax = OutputSyntax.SnpEff;

	@Parameter(names={"-R","--reference"},description= INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidxPath = null;
	

	/** a sequence with one or more altered base/amino-acid */
	private static class MutedSequence extends DelegateCharSequence
		{
		private final int pos;
		private final int refLen;
		private final String mut;
		MutedSequence(final CharSequence wild,final int pos,int refLen,final String mut)
			{
			super(wild);
			this.pos = pos;
			this.refLen = refLen;
			this.mut= mut;
			}
		@Override
		public int length() {
			return super.length() - this.refLen + this.mut.length();
			}
		
		@Override
		public char charAt(int i)
			{
			if(i < this.pos) return getDelegate().charAt(i);
			i-= this.pos;
			if(i< mut.length()) return this.mut.charAt(i);
			i-= this.mut.length();
			i+= this.refLen;//TODO check
			return getDelegate().charAt(this.pos+i);
			}		
		}
	
		
	
	private class Annotation
		{
		Transcript kg;
		final Allele alt;
		String exon_name="";
		String intron_name="";
		Integer position_cds=null;
		String wildAA="";
		String mutAA="";
		String wildCodon="";
		String mutCodon="";
		Integer position_protein=null;
		final Set<String> seqont=new HashSet<>();
		String transcriptId="";
		String geneId="";
		String geneName="";
		String bioType="";
		
		Annotation(final Allele alt) {
			this.alt = alt;
		}
		Annotation(final Allele alt,final Transcript transcript) {
			this(alt);
			this.kg = transcript;
			this.transcriptId = transcript.getId();
			this.geneId = transcript.getGene().getId();
			this.geneName = transcript.getGene().getGeneName();
			this.bioType = transcript.getGene().getGeneBiotype();
		}
		
		public String toString()
			{
			final StringBuilder b=new StringBuilder();
			switch(VCFPredictions.this.outputSyntax)
				{
				case Vep:
					{
					//Allele|Feature|Feature_type|Consequence|CDS_position|Protein_position|Amino_acids|Codons
					b.append(alt==null?"":alt.getBaseString());
					b.append('|');
					b.append(kg==null?"":kg.getId());
					b.append('|');
					b.append("Transcript");
					b.append('|');
					b.append(this.seqont.stream().map(T->T.replaceAll("[ ]", "_")).collect(Collectors.joining("&")));
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
					b.append(alt==null?"":alt.getBaseString());
					b.append('|');
					b.append(this.seqont.stream().map(T->T.replaceAll("[ ]", "_")).collect(Collectors.joining("&")));
					b.append('|');
					b.append("MODIFIER|");//impact 
					b.append(kg==null?"":kg.getGene().getGeneName());// gene name
					b.append('|');
					b.append(kg==null?"":kg.getGene().getId()); //gene id
					b.append('|');//
					b.append("transcript|");//Feature_Type
					b.append(kg==null?"":kg.getId());//Feature_ID
					b.append('|');
					b.append(kg==null?"":kg.getGene().getGeneBiotype());//Transcript_BioType
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
					b.append(kg==null?"":kg.getId());
					b.append('|');
					b.append(position_cds==null?"":String.valueOf(position_cds));
					b.append('|');
					b.append(position_protein==null?"":String.valueOf(position_protein));
					b.append('|');
					b.append(wildCodon.isEmpty() && mutCodon.isEmpty()?"":wildCodon+"/"+mutCodon);
					b.append('|');
					b.append(wildAA.isEmpty() && mutAA.isEmpty()?"":wildAA+"/"+mutAA);
					b.append("|");
					b.append(String.join("&",seqont));
					break;
					}
				}
			return b.toString();
			}
		
		}	
	
	private void loadGtf(final SAMSequenceDictionary dict) throws IOException
		{
		GtfReader in=null;
		try {
			this.transcriptTreeMap = new IntervalTreeMap<>();
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
			
			in = new GtfReader(this.gtfPath);
			in.setContigNameConverter(contigNameConverter);
			in.getAllGenes().
				stream().
				flatMap(G->G.getTranscripts().stream()).
				filter(G->G.hasCDS() && G.hasStrand()).
				forEach(g->{
				final int extend_gene_search = 5000; // because we want to set
														// SO:5KB_upstream_variant
				final Interval interval = new Interval(
						g.getContig(),
						Math.max(1, g.getTxStart() + 1 - extend_gene_search),
						g.getTxEnd() + extend_gene_search
						);
				List<Transcript> L= this.transcriptTreeMap.get(interval);
				if(L==null) {
					L=new ArrayList<>(2);
					this.transcriptTreeMap.put(interval, L);
				}
				L.add(g);
				});
			}
		finally {
			CloserUtil.close(in);
			}
		}
	private boolean isStop(char c)
		{
		return !Character.isLetter(c);
		}
	
	
	public static final String TAG="PRED";
	public static enum FORMAT1{TRANSCRIPT,CDSPOS,PROTPOS,CODON,AA,SEQONTOLOGY};
	
	private GenomicSequence getGenomicSequence(final String normalizedContig) {
		if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(normalizedContig))
			{
			this.transcriptId2cdna.clear();
			this.genomicSequence= new GenomicSequence(this.referenceGenome, normalizedContig);
			}
		return this.genomicSequence;
		}
	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator r, VariantContextWriter w) {
			try {
			this.referenceGenome= ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidxPath);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.referenceGenome);
			final ContigNameConverter contigNameConverter = ContigNameConverter.fromOneDictionary(dict);
			
			loadGtf(dict);
			
			final VCFHeader header= r.getHeader();
			final VCFHeader h2=new VCFHeader(header);
			JVarkitVersion.getInstance().addMetaData(this, h2);
			
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
	
	
			final RNASequenceFactory rnaSeqFactory = new RNASequenceFactory();
			rnaSeqFactory.setContigToGenomicSequence(S->getGenomicSequence(S));
			
			while(r.hasNext())
				{
				final VariantContext ctx= r.next();
				
				
				
				final String normalizedContig=contigNameConverter.apply(ctx.getContig());
				
				if(StringUtil.isBlank(normalizedContig)) {
					w.add(ctx);
					continue;
					}
				
				
				final List<Transcript> transcripts =  this.transcriptTreeMap.getOverlapping(new SimpleInterval(normalizedContig,ctx.getStart(),ctx.getEnd() )).
						stream().
						flatMap(L->L.stream()).
						collect(Collectors.toList());
				
				final List<Annotation> all_annotations = new ArrayList<>();
				final List<Allele> alternateAlleles;
				
			
				if(ctx.getNAlleles()<=1) {// not a variant, just REF
					alternateAlleles = Arrays.asList(Allele.NO_CALL);
					}
				else
					{
					alternateAlleles = ctx.getAlternateAlleles();
					}
				
				
				for(final Allele altAllele: alternateAlleles) {
					if(altAllele.isReference() || altAllele.equals(Allele.SPAN_DEL) || altAllele.equals(Allele.NON_REF_ALLELE) ) continue;
					
					/* intergenic ====================================================== */
					if(transcripts.isEmpty())
						{
						Transcript leftGene = null;
						String leftId="";
						String leftName="";
						
						for(Iterator<Transcript> iter=this.transcriptTreeMap.getOverlapping(new SimpleInterval(normalizedContig,1,ctx.getStart())).
								stream().
								flatMap(L->L.stream()).
								iterator();
								iter.hasNext();
								)
							{
							final Transcript t= iter.next();
							if(leftGene==null || leftGene.getEnd() < t.getEnd()) {
								leftGene = t;
								leftId= t.getGene().getId();
								leftName= t.getGene().getGeneName();
								}
							}
						
						Transcript rightGene = null;
						String rightId="";
						String rightName="";
	
						for(Iterator<Transcript> iter=this.transcriptTreeMap.getOverlapping(new SimpleInterval(normalizedContig,ctx.getEnd(),dict.getSequence(normalizedContig).getSequenceLength())).
								stream().flatMap(L->L.stream()).iterator();
								iter.hasNext();
								)
							{
							final Transcript t= iter.next();
							if(rightGene==null || t.getStart() < rightGene.getStart()) {
								rightGene = t;
								rightId= t.getGene().getId();
								rightName= t.getGene().getGeneName();
								}
							}
						
						//intergenic
						final Annotation annot = new Annotation(altAllele);
						annot.seqont.add("intergenic");
						annot.geneId=leftId+"-"+rightId;
						annot.geneName=leftName+"-"+rightName;
						all_annotations.add(annot);
						}
					else
						{
						
						//final int position1 =  ctx.getStart();
		
										
						for(final Transcript transcript:transcripts)
							{
							final Annotation annotation = new Annotation(altAllele,transcript);
							all_annotations.add(annotation);
							
							if(!transcript.overlaps(ctx)) {
								if( ((ctx.getEnd() < transcript.getStart() && transcript.isNegativeStrand()) ||
									  (ctx.getStart() > transcript.getEnd() && transcript.isPositiveStrand()))
									) {
									if(ctx.withinDistanceOf(transcript, 500)) {
										annotation.seqont.add("500B_downstream_variant");
										}
									else if(ctx.withinDistanceOf(transcript, 2_000))
										{
										annotation.seqont.add("2KB_downstream_variant");
										}
									}
								else if( ((ctx.getEnd() < transcript.getStart() && transcript.isPositiveStrand()) ||
										  (ctx.getStart() > transcript.getEnd() && transcript.isNegativeStrand())) 
									) {
									if(ctx.withinDistanceOf(transcript, 2_000)) {
										annotation.seqont.add("2KB_upstream_variant");
										}
									else if(ctx.withinDistanceOf(transcript, 5_000))
										{
										annotation.seqont.add("5KB_upstream_variant");
										}
									}
								continue;
								}
							
							if(CoordMath.encloses(ctx.getStart(), ctx.getEnd(), transcript.getStart(), transcript.getEnd())) {
								annotation.seqont.add("transcript_ablation");//TODO can be inversion ,etc...
								continue;
							}
							
							for(int side=0;side<2;++side) {
								final Optional<UTR> opt_utr= (side==0?transcript.getTranscriptUTR5():transcript.getTranscriptUTR3());
								if(!opt_utr.isPresent()) continue;
								final UTR utr = opt_utr.get();
								if(CoordMath.overlaps(utr.getStart(), utr.getEnd(), ctx.getStart(), ctx.getEnd()))
									{
									annotation.seqont.add(side==0?"5_prime_UTR_variant":"3_prime_UTR_variant");
									}
								
	 							}
							
							
							for(int side=0;side<2;++side) {
								final Optional<? extends ExonOrIntron> opt_ex ;
								if( side==0 ) {
									opt_ex = transcript.getExons().stream().filter(E->E.overlaps(ctx)).findFirst();
									}
								else
									{
									opt_ex = transcript.getIntrons().stream().filter(E->E.overlaps(ctx)).findFirst();
									}
								if(!opt_ex.isPresent()) continue;
								final ExonOrIntron ei = opt_ex.get();
								if(side==0) {
									if(transcript.isNonCoding()) annotation.seqont.add("non_coding_transcript_exon_variant");
									}
								else
									{
									if(transcript.isNonCoding()) annotation.seqont.add("non_coding_transcript_intron_variant");
									annotation.seqont.add("intron");
									}
								
								if(ctx.getStart()==ctx.getEnd() && ei.isSplicing(ctx.getStart()))
									{
									if(ei.isSplicingAcceptor(ctx.getStart()))
										{
										annotation.seqont.add("splice_acceptor"); //SPLICING_ACCEPTOR
										}
									else  if(ei.isSplicingDonor(ctx.getStart()))
										{
										annotation.seqont.add("splice_donor"); // SPLICING_DONOR
										}
									else //??
										{
										annotation.seqont.add("splicing_variant");
										}
									}
								}
		
							
							final StructuralVariantType svType = ctx.getStructuralVariantType();
							if(svType!=null) {
								
								
								
								continue;
							}
							
							if(transcript.isNonCoding()) {
								//TODO
								annotation.seqont.add("non_coding_transcript_variant");
								continue;
							}
							
							RNASequence cDNA = this.transcriptId2cdna.get(transcript.getId());
							if(cDNA==null) {
								cDNA = rnaSeqFactory.getCodingRNA(transcript);
								 this.transcriptId2cdna.put(transcript.getId(), cDNA);
								}
							final OptionalInt opt_pos_cdna0 = cDNA.convertGenomic0ToRnaIndex0(ctx.getStart()-1);
							if(!opt_pos_cdna0.isPresent()) continue;
							final int pos_cdna0 = opt_pos_cdna0.getAsInt();
		            		final int pos_aa = pos_cdna0/3;
		
							
							final GeneticCode geneticCode=GeneticCode.getStandard();
		            							
							if(AcidNucleics.isATGC(altAllele.getBaseString()))
								{
								
		 						String bases = altAllele.getBaseString().toUpperCase();
								if(transcript.isNegativeStrand()) {
									bases = AcidNucleics.reverseComplement(bases);
									}
								final MutedSequence mutRNA = new MutedSequence(cDNA,pos_cdna0,ctx.getReference().length(),bases);
								final PeptideSequence<CharSequence> wildProt = PeptideSequence.of(cDNA,geneticCode);
								final PeptideSequence<CharSequence> mutProt = PeptideSequence.of(mutRNA,geneticCode);
								
			            		final int mod= pos_cdna0%3;
			            		annotation.wildCodon=(""+
			            			cDNA.charAt(pos_cdna0-mod+0)+
			            			cDNA.charAt(pos_cdna0-mod+1)+
			            			cDNA.charAt(pos_cdna0-mod+2)
			            			);
			            		annotation.mutCodon=(""+
			            			mutRNA.charAt(pos_cdna0-mod+0)+
			            			mutRNA.charAt(pos_cdna0-mod+1)+
			            			mutRNA.charAt(pos_cdna0-mod+2)
			            			);
			            		annotation.position_protein=(pos_aa+1);
			            		annotation.wildAA=String.valueOf(wildProt.charAt(pos_aa));
			            		annotation.mutAA=(String.valueOf(mutProt.charAt(pos_aa)));
			            		
			            		
				    			if(isStop(wildProt.charAt(pos_aa)) &&
				    			   !isStop(mutProt.charAt(pos_aa)))
				    				{
				    				annotation.seqont.add("stop_lost");
				    				}
				    			else if( !isStop(wildProt.charAt(pos_aa)) &&
				    				 isStop(mutProt.charAt(pos_aa)))
				    				{
				    				annotation.seqont.add("stop_gained");
				    				}
				    			else if(wildProt.charAt(pos_aa)==mutProt.charAt(pos_aa))
				    				{
				    				annotation.seqont.add("synonymous");
				    				}
				    			else
				    				{
				    				annotation.seqont.add("missense_variant");
				    				}								
								}
							}
						}
					}
				
			
				
				final Set<String> info=new HashSet<String>(all_annotations.size());
				for(final Annotation a:all_annotations)
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
			return 0;
			}
		catch(final Throwable err ) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(w);
			CloserUtil.close(r);
			CloserUtil.close(this.referenceGenome);
			}
		}
	
	
	public static void main(final String[] args)
		{
		new VCFPredictions().instanceMainWithExit(args);
		}
	}
