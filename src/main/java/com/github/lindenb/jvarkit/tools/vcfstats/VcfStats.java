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
package com.github.lindenb.jvarkit.tools.vcfstats;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.RangeOfIntegers;
import com.github.lindenb.jvarkit.tools.burden.MafCalculator;
import com.github.lindenb.jvarkit.tools.lumpysv.LumpyConstants;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
/*
BEGIN_DOC

## Tip: Adding a new key in the INFO field

Using vcffilterjs :



the script:

```
var ArrayList = Java.type("java.util.ArrayList");
var VariantContextBuilder = Java.type("htsjdk.variant.variantcontext.VariantContextBuilder");


function addInfo(v)
	{
	var vcb = new VariantContextBuilder(v);
	var atts = new ArrayList();
	atts.add(v.getType().name()+ (variant.isFiltered()?"_FILTERED":"_UNFILTERED"));
	atts.add(v.getType().name()+ (variant.hasID()?"_ID":"_NOID"));
	vcb.attribute("MYKEY",atts);
	return vcb.make();
	}


addInfo(variant);
```

run the program, but first use awk to insert the new INFO definition for 'MYKEY'

```
cat input.vcf |\
	awk '/^#CHROM/ {printf("##INFO=<ID=MYKEY,Number=.,Type=String,Description=\"My key\">\n");} {print;}' |\
	java -jar dist/vcffilterjs.jar -f script.js 
```


## Example

```
$ java -jar $< -o tmp --soterms SO:0001818 --soterms SO:0001819 -K ucsc/hg19/database/knownGene_noPrefix.txt.gz input.vcf
$ (cd tmp && make) 
$ ls tmp/
ALL.affectedSamples.png     ALL.countDepthBySample.tsv      ALL.countDistances.png  ALL.geneLoc.tsv  ALL.predictionsBySample.png  ALL.sample2gtype.tsv  Makefile
ALL.affectedSamples.tsv     ALL.countDepth.png              ALL.countDistances.tsv  ALL.maf.png      ALL.predictionsBySample.tsv  ALL.transvers.png
ALL.countAltAlleles.png     ALL.countDepth.tsv              ALL.countIndelSize.png  ALL.maf.tsv      ALL.predictions.png          ALL.transvers.tsv
ALL.countAltAlleles.tsv     ALL.countDistancesBySample.png  ALL.countIndelSize.tsv  ALL.mendel.png   ALL.predictions.tsv          ALL.variant2type.png
ALL.countDepthBySample.png  ALL.countDistancesBySample.tsv  ALL.geneLoc.png         ALL.mendel.tsv   ALL.sample2gtype.png         ALL.variant2type.tsv

$ head -n3 tmp/*.tsv
==> tmp/ALL.affectedSamples.tsv <==
1/57	182265
2/57	98512
3/57	67449

==> tmp/ALL.countAltAlleles.tsv <==
2	11699
3	2406
4	679

==> tmp/ALL.countDepthBySample.tsv <==
Sample	[-Inf/0[	[0/10[	[10/20[	[20/30[	[30/50[	[50/100[	[100/200[	[200/Inf[
10_T1245	0	305229	97739	61972	75791	72779	14008	686
11AG09	0	400147	62261	39411	55254	86883	86192	46749

==> tmp/ALL.countDepth.tsv <==
[0/10[	19435
[10/20[	95369
[20/30[	92378

==> tmp/ALL.countDistancesBySample.tsv <==
Sample	[-Inf/0[	0	1	2	3	4	5	6	7	8	9	[10/20[	[20/100[	[100/200[	[200/300[	[300/400[	[400/500[	[500/1000[	[1000/Inf[
S1	0	0	2643	1755	1491	1430	1156	1106	893	858	867	6826	26696	18252	12411	9348	7270	20975	85769
S2	0	0	3632	2397	1967	1855	1516	1425	1246	1179	1109	8917	36006	25455	18652	14951	11912	37828	103341

==> tmp/ALL.countDistances.tsv <==
1	16275
2	11184
3	8505

==> tmp/ALL.countIndelSize.tsv <==
2	59420
3	22727
4	10216

==> tmp/ALL.geneLoc.tsv <==
first_exon	51180
internal_exon	70074
last_exon	75341

==> tmp/ALL.maf.tsv <==
0.029411764705882353	0.027777777777777776
0.5	0.43333333333333335
0.05263157894736842	0.0

==> tmp/ALL.mendel.tsv <==
Sample	synonymous_variant	protein_altering_variant
10_T1245	0	0
11AG09	170	298

==> tmp/ALL.predictionsBySample.tsv <==
Sample	synonymous_variant	protein_altering_variant
S1	13453	14527
S2	12820	14077

==> tmp/ALL.predictions.tsv <==
protein_altering_variant	59516
synonymous_variant	47454

==> tmp/ALL.sample2gtype.tsv <==
Sample	NO_CALL	HOM_REF	HET	HOM_VAR	UNAVAILABLE	MIXED
S1	345776	429002	98825	100945	0	0
S2	196925	504211	132576	140836	0	0

==> tmp/ALL.transvers.tsv <==
TYPE	TRANSITION	TRANSVERSION
ALL	581537	133472
CDS	533340	120720

==> tmp/ALL.variant2type.tsv <==
Type	Count
INDEL	126219
SNP	846677

```


END_DOC
 */
@Program(name="vcfstats",
	description="Produce VCF statitics",
	keywords={"vcf","stats","burden","gnuplot"}
	)
public class VcfStats extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfStats.class).make();

	@Parameter(names={"-o","--output"},description="output Directory or zip file. The output contains the data files as well as a Makefile to convert the data files to graphics using gnuplot.",required=true)
	private File outputFile = null;
	
	@Parameter(names={"-K","-kg","--knownGenes"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String kgFile = null;
	private IntervalTreeMap<List<KnownGene>> knownGeneTreeMap=null;
	
	@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile = null;
	private Pedigree pedigree = null;

	
	@Parameter(names={"--prefix"},description="File/zip prefix")
	private String prefix = "";
	
	@Parameter(names={"--vckey"},description="Variant Context Key. if defined, I will look at this key in the INFO column and produce a CASE/CTRL graf for each item. If undefined, I will produce a default graph with all variant")
    private String mafKeyInINFO = null;
	@Parameter(names={"-mafTag","--mafTag"},description="Do not calculate MAF for controls, but use this tag to get Controls' MAF")
	private String controlTag =null;
	@Parameter(names={"-nchr","--nocallishomref"},description="treat no call as HomRef")
	boolean no_call_is_homref=false;
	@Parameter(names={"-tee","--tee"},description="output the incoming vcf to stdout. Useful to get intermediary stats in a pipeline")
	boolean tee=false;
	@Parameter(names={"--trancheAffected"},description="tranches for the number of affected. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers affectedTranches = new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,20,50,100,200,300,400,500,1000);
	@Parameter(names={"--trancheDP"},description="tranches for the DEPTH. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers depthTranches =new RangeOfIntegers(0,10,20,30,50,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,10_000,20_000,30_000,40_000,50_000,100_000);
	@Parameter(names={"--trancheIndelSize"},description="tranches for the Indel size "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers indelTranches =new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,15,20,50,100);
	@Parameter(names={"--trancheAlts"},description="tranches for the number of ALTs. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers altTranches =new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10);
	@Parameter(names={"--trancheDistance"},description="tranches for the distance between the variants. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers distanceTranches =new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,20,100,200,300,400,500,1000);
	@Parameter(names={"-so","--soterms"},description="Sequence ontology Accession to observe. VCF must be annotated with SNPEFF or VEP. e.g: \"SO:0001818\" (protein altering variant) \"SO:0001819\" (synonymouse variant)")
	private Set<String> sequenceOntologyTermsStr=new HashSet<>();
	@Parameter(names={"--disableMAFPlot"},description="Disable MAF plot")
	private boolean disableMAFPlot=false;
	@Parameter(names={"--disableGTConcordance"},description="Disable Plot Sample vs Sample Genotypes (Faster...)")
	private boolean disableGenotypeConcordance=false;
	@Parameter(names={"--binSize"},description="[20170718] When plotting data over a genome, divide it into 'N' bp.")
	private int binSize = 1_000_000;
	
	private ArchiveFactory archiveFactory=null;
	/** the SAMSequenceDictionary used to sort reference */
	private SAMSequenceDictionary the_dictionary = null;
	/** list of samples in order*/
	private List<String> sampleNamesInOrder = Collections.emptyList();
	
	private final Function<String, Integer> contig2tid = (S)->{
		final int tid = the_dictionary.getSequenceIndex(S);
		if(tid<0) throw new JvarkitException.ContigNotFoundInDictionary(S, the_dictionary);
		return tid;
		};
	
	private final Comparator<String> contigComparator = (S1,S2) -> {
		if(S1.equals(S2)) return 0;
		if(the_dictionary==null) {
			return S1.compareTo(S2);
			} else {
				return contig2tid.apply(S1) - contig2tid.apply(S2);
			}
		};
	
	private class ContigBin	implements Locatable,Comparable<ContigBin> {
		final String contig;
		final int pos;
		ContigBin(final String contig,final int pos) {
			this.contig = contig;
			this.pos = pos - pos%VcfStats.this.binSize;
		}
		
		@Override
		public String getContig() {
			return this.contig;
			}
		@Override
		public int getStart() {
			return pos+1;
			}
		@Override
		public int getEnd() {
			return pos + VcfStats.this.binSize;
			}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + contig.hashCode();
			result = prime * result + pos;
			return result;
			}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null || getClass() != obj.getClass())
				return false;
			final ContigBin other = (ContigBin) obj;
			return this.pos == other.pos && 
				  this.contig.equals(other.contig) ; 
			}
		@Override
		public int compareTo(final ContigBin o) {
			int i = VcfStats.this.contigComparator.compare(this.getContig(), o.getContig());
			if(i!=0) return i;
			return pos - o.pos;
			}
		@Override
		public String toString() {
			return getContig()+";"+this.getStart()+"-"+this.getEnd();
			}
		}
	
	//define a pair of samples
	private class SamplePair
		{
		final int sample1;
		final int sample2;
		
		SamplePair(final int sample1,final int sample2)
			{
			if(sample1 < sample2)
				{
				this.sample1 = sample1;
				this.sample2 = sample2;
				}
			else
				{
				this.sample1 = sample2;
				this.sample2 = sample1;
				}
			}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + this.sample1;
			result = prime * result + this.sample2;
			return result;
			}
		@Override
		public boolean equals(final Object obj) {
			if (this == obj)  return true;
			if (obj == null)  return false;
			if (!(obj instanceof SamplePair)) return false;
			final SamplePair other = (SamplePair) obj;
			return other.sample1==this.sample1 && 
				   other.sample2==this.sample2; 
			}
		@Override
		public String toString() {
			return  VcfStats.this.sampleNamesInOrder.get(sample1)+"/"+
					VcfStats.this.sampleNamesInOrder.get(sample2);
			}
		}
	
	
	private class PlotMaf implements Closeable
		{
		//final String title;
		final String filename;
		final PrintWriter pw;
		PlotMaf(final String title)
			{
			//this.title=title;
			this.filename = VcfStats.this.prefix + title+".maf.tsv";
			try {
				this.pw = archiveFactory.openWriter(this.filename);
			} catch (IOException e) {
				throw new RuntimeIOException(e);
			}
			}
		void plot(double xcas,double yctrl) {
			this.pw.print(xcas);
			this.pw.print('\t');
			this.pw.print(yctrl);
			this.pw.print('\n');
			}	
		@Override
		public void close() throws IOException {
			this.pw.flush();
			this.pw.close();
			}
		}
	

		private final Function<VariantContext, Set<String>> variantToCategoryKeys = VC ->{
			final Set<String> set = new HashSet<>();
			set.add("ALL");
			if( !(mafKeyInINFO==null || mafKeyInINFO.trim().isEmpty())) {
				set.addAll( VC.getAttributeAsList(mafKeyInINFO).
					stream().
					filter(O->!(O==null || ".".equals(O))).
					map(O->String.valueOf(O)).
					collect(Collectors.toSet()));
				}
			return set;
			};

		
	// where are the mutations ordered by priority 
	private enum GeneLocation
		{
		first_exon,
		internal_exon,
		last_exon,
		not_in_gene
		}
		
	/** statistics common to variant and samples */
	private abstract class AbstractStat
		{
		final Counter<RangeOfIntegers.Range> countDepth = new Counter<>();
		final Counter<RangeOfIntegers.Range> countDistances = new Counter<>();
		final Counter<NucleicAcidChange> nucleicAcidChanges = new Counter<>();
		final Counter<GeneLocation> geneLocations = new Counter<>();
		final Counter<String> consequences = new Counter<>();
		final Counter<String> variantsPerContigs = new Counter<>();
		final Counter<ContigBin> countBins = new Counter<>();
		final Counter<StructuralVariantType> countStructuralVariations = new Counter<>();
		
		protected ContigPosRef prevCtx=null;

		protected void visitForDistance(final VariantContext ctx)
			{
			//distance
			final ContigPosRef contigPosRef = new ContigPosRef(ctx);
			if(prevCtx!=null && prevCtx.getContig().equals(contigPosRef.getContig()) && prevCtx.getStart() <= contigPosRef.getStart())
				{
				final int distance = contigPosRef.getStart() - this.prevCtx.getStart();
				this.countDistances.incr(VcfStats.this.distanceTranches.getRange(distance));
				}
			prevCtx=contigPosRef;
			}
		
		protected void visitForConsequences(final VariantContext ctx,Collection<SequenceOntologyTree.Term> terms,final VcfTools vcfTools)
			{
			for(final SequenceOntologyTree.Term t:terms)
				{
				if(vcfTools.hasSequenceOntologyTerm(ctx, t))
					{
					this.consequences.incr(t.getLabel());
					}	
				}
			}

		
		protected void visitForGeneLocation(final VariantContext ctx,final List<KnownGene> kgs)
			{
			GeneLocation loc= GeneLocation.not_in_gene;
			for(final KnownGene kg:kgs)
				{
				for(final KnownGene.Exon exon:kg.getExons())
					{
					if(!exon.contains(ctx.getStart()-1)) continue;
					GeneLocation loc2;
					if(exon.isFirstExon())
						{
						loc = GeneLocation.first_exon;//best 
						break;
						}
					else if(exon.isLastExon())
						{
						 loc2 = GeneLocation.last_exon;
						}
					else
						{
						loc2 = GeneLocation.internal_exon;
						}
						
					if(loc2.compareTo(loc)<0) loc=loc2;
					}
				}
			if(loc!=GeneLocation.not_in_gene)
				{
				this.geneLocations.incr(loc);
				}
			}
		
		}
	
	
	private enum NucleicAcidChange
		{
		transition,
		transversion,
		transition_in_cds,
		transversion_in_cds
		}
	
	
	private class VariantStats extends AbstractStat
		{
		private final VcfTools vcfTools;
		private final Set<SequenceOntologyTree.Term> sequenceOntologyTermsToObserve=new HashSet<>();
		private final String key;
		private PlotMaf mafPlotter= null;
		private final Set<String> affectedSamples;
		private final Set<String> unaffectedSamples;
		private final Map<String,SampleStat> sample2stats = new TreeMap<>();
		final Counter<VariantContext.Type> countTypes = new Counter<>();
		final Counter<RangeOfIntegers.Range> countAffectedSamples = new Counter<>();
		final Counter<RangeOfIntegers.Range> countAltAlleles = new Counter<>();
		final Counter<RangeOfIntegers.Range> countIndelSize = new Counter<>();
		final Counter<SamplePair>  genotypeConcordance = new Counter<>();

		private int countVariants=0;
		
		/** stats for Samples */
		private class SampleStat extends AbstractStat
			{
			final Counter<GenotypeType> countTypes = new Counter<>();
			final Counter<String> countMendelianViolations = new Counter<>();
			final String sampleName;
			final Pedigree.Person pedireePerson;
			SampleStat(final String sampleName) {
				this.sampleName = sampleName;
				this.pedireePerson= VcfStats.this.pedigree.getPersonById(sampleName);
				}
			
			public void visit(final VariantContext ctx,final List<KnownGene> knownGenes) {
				final Genotype genotype = ctx.getGenotype(this.sampleName);
				if(genotype==null) return;
				this.countTypes.incr(genotype.getType());
				if(ctx.isVariant() && genotype.isCalled() && !genotype.isHomRef())
					{
					this.variantsPerContigs.incr(ctx.getContig());
					this.countBins.incr(new ContigBin(ctx.getContig(),ctx.getStart()));
					}
				
				if( this.pedireePerson !=null && 
					this.pedireePerson.hasAtLeastOneParent() &&
					VariantStats.this.vcfTools.isMendelianIncompatibility(ctx, this.pedireePerson))
					{
					for(final SequenceOntologyTree.Term t:VariantStats.this.sequenceOntologyTermsToObserve)
						{
						if(vcfTools.hasSequenceOntologyTerm(ctx, t))
							{
							this.countMendelianViolations.incr(t.getLabel());
							}
						}
					}
				final StructuralVariantType structuralVariantType= ctx.getStructuralVariantType();
				if(structuralVariantType!=null)
					{
					if( 
						(genotype.isCalled() && !genotype.isHomRef()) ||
						(LumpyConstants.isLumpyVariant(ctx) && genotype.hasExtendedAttribute("SU") && genotype.getAttributeAsInt("SU", 0)>0)
						)
						{
						this.countStructuralVariations.incr(structuralVariantType);
						}
					
					}
				
				
				if(genotype.hasDP())
					{
					final int dp = genotype.getDP();
					if(dp>=0)
						{
						this.countDepth.incr(VcfStats.this.depthTranches.getRange(dp));
						}
					}
				
				if(genotype.isHomVar() || genotype.isHet())
					{
					visitForDistance(ctx);
					visitForGeneLocation(ctx,knownGenes);
					visitForConsequences(ctx,
							VariantStats.this.sequenceOntologyTermsToObserve,
							VariantStats.this.vcfTools
							);

					}
				}
			public void finish(final PrintWriter makefileWriter) throws IOException
				{
				
				}
	
			}

		
		
		
		VariantStats(final String key,final VCFHeader header) {
			this.key = key;
			this.vcfTools = new VcfTools(header);
			
			
			this.sequenceOntologyTermsToObserve.addAll(
				VcfStats.this.sequenceOntologyTermsStr.stream().
			 	filter(S->!S.trim().isEmpty()).
			 	map(S->SequenceOntologyTree.createDefault().getTermByAcn(S)).
			 	collect(Collectors.toSet())
			 	);
			
			
			
			for(final String sn: VcfStats.this.sampleNamesInOrder)
				{
				this.sample2stats.put(sn,new SampleStat(sn));
				}
			
			this.affectedSamples = 
					VcfStats.this.pedigree.getPersons().stream().
						filter(P->P.isAffected()).
						map(P->P.getId()).
						filter(S->VcfStats.this.sampleNamesInOrder.contains(S)).
						collect(Collectors.toSet())
						;
			
			this.unaffectedSamples = 
					VcfStats.this.pedigree.getPersons().stream().
						filter(P->P.isUnaffected()).
						map(P->P.getId()).
						filter(S->VcfStats.this.sampleNamesInOrder.contains(S)).
						collect(Collectors.toSet())
						;
			// genotype concordance
			if(!VcfStats.this.disableGenotypeConcordance) {
				final int n_samples = VcfStats.this.sampleNamesInOrder.size();
				for(int x=0;x +1 < n_samples;++x)
					{
					for(int y= x +1; y < n_samples;++y)
						{
						this.genotypeConcordance.initializeIfNotExists(new SamplePair(x,y));
						}
					}
				}
			}
		
		public void visit(final VariantContext ctx) {
			this.countVariants++;
			this.countTypes.incr(ctx.getType());
			if(ctx.isVariant())
				{
				this.variantsPerContigs.incr(ctx.getContig());
				this.countBins.incr(new ContigBin(ctx.getContig(),ctx.getStart()));
				}
			
			
			final List<KnownGene> knownGenes =  VcfStats.this.getOverlappingKnownGenes(ctx);			

			visitForGeneLocation(ctx,knownGenes);
			for(final SampleStat st: this.sample2stats.values()) st.visit(ctx,knownGenes);
			
			//distance
			visitForDistance(ctx);
			
			final StructuralVariantType structuralVariantType= ctx.getStructuralVariantType();
			if(structuralVariantType!=null)
				{
				this.countStructuralVariations.incr(structuralVariantType);
				}

			
			final List<Allele> alternates = ctx.getAlternateAlleles();

			
			/** consequences */
			visitForConsequences(ctx,this.sequenceOntologyTermsToObserve,this.vcfTools);
			
			/** transvertion / transition */
			if(alternates.size()==1 )
				{
				boolean in_cds=false;
				if(VcfStats.this.knownGeneTreeMap!=null)
					{
					in_cds = knownGenes.stream().
							filter(K->!((K.getTxStart()+1) > ctx.getEnd() || (K.getTxEnd()) < ctx.getStart()  )).
							flatMap(K->K.getExons().stream()).
							filter(E->!((E.getStart()+1 > ctx.getEnd() ))).
							findAny().
							isPresent()
							;
					}
				
				final Character refChar=asSimpleATGC(ctx.getReference());
				final Character altChar=asSimpleATGC(alternates.get(0));
				if(isTransition(refChar,altChar)) {
					this.nucleicAcidChanges.incr(NucleicAcidChange.transition);
					if(VcfStats.this.knownGeneTreeMap!=null && in_cds) {
						this.nucleicAcidChanges.incr(NucleicAcidChange.transition_in_cds);
						}
					}	
				else if(isTransversion(refChar,altChar)) {
					this.nucleicAcidChanges.incr(NucleicAcidChange.transversion);
					if(VcfStats.this.knownGeneTreeMap!=null && in_cds) {
						this.nucleicAcidChanges.incr(NucleicAcidChange.transversion_in_cds);
						}
					}	
				}
			
			/* position in the genes */
			
			
			/* can we compute the MAF ? we need (non)affected samples*/
			if(!VcfStats.this.disableMAFPlot && !(alternates.isEmpty() || this.unaffectedSamples.isEmpty() || this.affectedSamples.isEmpty()))
				{
				for(int alt_idx=0;alt_idx < alternates.size();++alt_idx) {
					final Allele alt = alternates.get(alt_idx);
					final Double mafs[]={null,null};
					
					for(int i=0;i< 2;++i)
						{
						if(i==1 && VcfStats.this.controlTag!=null)
							{
							if(ctx.hasAttribute(VcfStats.this.controlTag)) {
								try 
									{
									final List<Double> dvals =ctx.getAttributeAsDoubleList(VcfStats.this.controlTag, Double.NaN);
									if(alt_idx< dvals.size() && dvals.get(alt_idx)!=null) {	
									    final double d= dvals.get(alt_idx);
										if(!Double.isNaN(d) && d>=0 && d<=1.0) mafs[1]=d;
										}
									}
								catch(NumberFormatException err)
									{
									}
								}
							}
						else
							{
							final MafCalculator mafCalculator = new MafCalculator(alt, ctx.getContig());
							mafCalculator.setNoCallIsHomRef(no_call_is_homref);
							for(Pedigree.Person person: (i==0?pedigree.getAffected():pedigree.getUnaffected()))
								{
								final Genotype genotype = ctx.getGenotype(person.getId());
								if(genotype==null) continue;
								mafCalculator.add(genotype, person.isMale());
								}
							if(!mafCalculator.isEmpty())
								{
								mafs[i]=mafCalculator.getMaf();
								}
							}
						}
					if(mafs[0]==null || mafs[1]==null) continue;
					
					
						
						if(this.mafPlotter==null)
							{
							//it's a new plotter
							this.mafPlotter = new PlotMaf(key);
							
							//add makefile stuff
							/**
							final String png= "$(patsubst %.tsv,%.png,"+plotter.filename+")";
							makefileWriter.println("ALL_TARGETS+=" + png);
							makefileWriter.println(png+":"+plotter.filename+" "+generic_maf_gnuplot_filename);
							makefileWriter.println("\tsed -e '%__OUTPUT__%$@%g' -e '%__INPUT__%$<%g'  -e '%__TITLE__%"+mafkey+"%g' $(word 2,$^) | gnuplot");
							*/
							}
						this.mafPlotter.plot(mafs[0], mafs[1]);
						
					
					}//end of loop over ALT
				} // end of MAF
			
			
			this.countAffectedSamples.incr(
					VcfStats.this.affectedTranches.getRange(
						(int)ctx.getGenotypes().stream().
							filter(G->G.isCalled() && !(G.isHomRef() || G.isFiltered() )).
							count()	)
					);
			if(ctx.hasAttribute(VCFConstants.DEPTH_KEY))
				{
				int dp = ctx.getAttributeAsInt(VCFConstants.DEPTH_KEY, -1);
				if(dp>=0)
					{
					this.countDepth.incr(VcfStats.this.depthTranches.getRange(dp));
					}
				}
			if(ctx.isIndel())
				{
				final int longest = ctx.getAlleles().stream().
						filter(A->!(A.isSymbolic() || A.equals(Allele.SPAN_DEL))).
						mapToInt(A->A.length()).max().orElse(0);
				this.countIndelSize.incr(VcfStats.this.indelTranches.getRange(longest));
				}	
			this.countAltAlleles.incr(VcfStats.this.altTranches.getRange(alternates.size()));
			
			// genotype concordance

			if(!VcfStats.this.disableGenotypeConcordance) {
				for(int x=0;x < ctx.getNSamples();++x)
					{
					final Genotype g1 = ctx.getGenotype(x);
					if(!g1.isCalled()) continue;
					for(int y= x ; y < ctx.getNSamples();++y)
						{
						final Genotype g2 = ctx.getGenotype(y);
						if(!g2.isCalled()) continue;
						if(g1.sameGenotype(g2))
							{
							this.genotypeConcordance.incr(new SamplePair(x,y));
							}
						}
					}
				}
			
			}
		private String toTsv(final String filename)
			{
			return VcfStats.this.prefix+this.key+"."+filename+".tsv";
			}

		
		private String toPng(final String filename)
			{
			if(!filename.endsWith(".tsv")) throw new IllegalArgumentException(filename);
			return "$(patsubst %.tsv,%.png,"+filename+")";
			}
		
		/** output results */
		public void finish(final PrintWriter makefileWriter) throws IOException
		{

		if(!this.countTypes.isEmpty())
			{
			final String filename = toTsv("variant2type");

			PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
			pw.println("Type\tCount");

			for(final VariantContext.Type type: this.countTypes.keySet())
				{
				pw.println(type.name()+"\t"+this.countTypes.count(type));
				}
			pw.flush();
			pw.close();
			final String png= toPng(filename);
			makefileWriter.println("ALL_TARGETS+=" + png);
			makefileWriter.println(png+":"+filename);
			makefileWriter.println("\techo 'set key autotitle columnheader;"
					+ "set ylabel \"Count\";"
					+ "set yrange [0:];"
					+ "set size  ratio 0.618;"
					+ "set title \"Variant Types\";"
					+ "set style fill solid border -1;"
					+ "set key  off;set datafile separator \"\t\";"
					+ "set auto x;"
					+ "set xtic rotate by 90 right;"
					+ "set style histogram;"
					+ "set style data histogram;"
					+ "set terminal png truecolor;"
					+ "set output \"$@\";"
					+ "plot \"$<\" using 2:xtic(1) ti \"Variant Type\";' | "
					+ "gnuplot");
			}

		if(!this.countBins.isEmpty()) {
			final String filename = toTsv("genscan");
			PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
			final String png= toPng(filename);
			makefileWriter.println("ALL_TARGETS+=" + png);
			makefileWriter.println(png+":"+filename);

			if(VcfStats.this.the_dictionary != null ) {
				long gpos = 0L;
				for(final SAMSequenceRecord ssr: VcfStats.this.the_dictionary.getSequences()) {
					int x = 0;
					while( x < ssr.getSequenceLength() ) {
						final long c = this.countBins.count(new ContigBin(ssr.getSequenceName(), x));
						if(c>0L) {
							pw.print(String.valueOf(gpos+x));
							pw.print('\t');
							pw.print(c);
							pw.print('\t');
							pw.print(ssr.getSequenceIndex());
							pw.print('\t');
							pw.print(ssr.getSequenceName()+":"+(x));
							pw.println();
							}
						x+= VcfStats.this.binSize;
						}
					gpos += ssr.getSequenceLength();
					}
				
					makefileWriter.println("\techo '"
						+ "set ylabel \"Count\";"
						+ "set xlabel \"Position Genome\";"
						+ "set yrange [0:];"
						+ "set title \"Variants over the Genome\";"
						+ "set style fill solid border -1;"
						+ "set key  off;"
						+ "set datafile separator \"\t\";"
						+ "set auto x;"
						+ "set xtic rotate by 90 right;"
						+ "set terminal png truecolor;"
						+ "set output \"$@\";"
						+ "plot \"$<\" using 1:2:3 w points linecolor variable' | "
						+ "gnuplot");
				
					} 
				else
					{
					for(final ContigBin ctgBin : new TreeSet<>(this.countBins.keySet())) {
						pw.print(ctgBin.toString());
						pw.print('\t');
						pw.print(this.countBins.count(ctgBin));
						pw.println();
						}
					
					makefileWriter.println("\techo '"
							+ "set ylabel \"Count\";"
							+ "set xlabel \"Position Genome\";"
							+ "set yrange [0:];"
							+ "set title \"Variants over the Genome\";"
							+ "set style fill solid border -1;"
							+ "set key  off;"
							+ "set datafile separator \"\t\";"
							+ "set auto x;"
							+ "set xtic rotate by 90 right;"
							+ "set style histogram;"
							+ "set style data histogram;"
							+ "set terminal png truecolor;"
							+ "set output \"$@\";"
							+ "plot \"$<\" using 2:xtic(1) ti \"Position\";' | "
							+ "gnuplot");
					}
			pw.flush();
			pw.close();
			}
		
		if(!this.variantsPerContigs.isEmpty())
			{
			final String filename = toTsv("variant2contigs");

			PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
			pw.println("Contig\tCount");

			for(final String contig: this.variantsPerContigs.keySetDecreasing())
				{
				pw.println(contig+"\t"+this.variantsPerContigs.count(contig));
				}
			pw.flush();
			pw.close();
			final String png= toPng(filename);
			makefileWriter.println("ALL_TARGETS+=" + png);
			makefileWriter.println(png+":"+filename);
			makefileWriter.println("\techo 'set key autotitle columnheader;"
					+ "set ylabel \"Count\";"
					+ "set yrange [0:];"
					+ "set size  ratio 0.618;"
					+ "set title \"Variant per Contigs (N="+this.variantsPerContigs.getTotal()+")\";"
					+ "set style fill solid border -1;"
					+ "set key  off; set datafile separator \"\t\";"
					+ "set auto x;"
					+ "set xtic rotate by 90 right;"
					+ "set style histogram;"
					+ "set style data histogram;"
					+ "set terminal png truecolor;"
					+ "set output \"$@\";"
					+ "plot \"$<\" using 2:xtic(1) ti \"Contig\";' | "
					+ "gnuplot");
			}

		if(!this.sample2stats.isEmpty())
			{
			final String filename=toTsv("sample2contig");
			final Set<String> contigs = new TreeSet<String>( VcfStats.this.contigComparator);
			contigs.addAll(this.sample2stats.values().
					stream().
					flatMap(S->S.variantsPerContigs.keySet().stream()).
					collect(Collectors.toSet())
					);
			
			PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
			pw.println("Sample\t"+ String.join("\t",contigs));
			for(final String sample: this.sample2stats.keySet())
				{
				pw.print(sample);
				for(final String contig: contigs)
					{
					pw.print("\t"+this.sample2stats.get(sample).variantsPerContigs.count(contig));
					}
				pw.println();
				}
			pw.flush();
			pw.close();
			
			final String png=toPng(filename);
			makefileWriter.println("ALL_TARGETS+=" + png);
			makefileWriter.println(png+":"+filename);

			makefileWriter.print("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT} ;"
					+ "set title \"Sample/Contig\";"
					+ "set xlabel \"Sample\";"
					+ "set key autotitle columnhead;"
					+ "set xtic rotate by 90 right;"
					+ "set ylabel \"Count\";"
					+ "set yrange [0:];"
					+ "set key invert reverse Left outside;"
					+ "set datafile separator \"\t\";"
					+ "set style fill solid border -1;"
					+ "set style data histograms;set style histogram rowstacked;"
					+ "set boxwidth 0.95;set output \"$@\";"
					+ "plot \"$<\" using 2:xtic(1)");
			int k=2;
			for(final String contig: contigs)
				{
				makefileWriter.print((k==2?"":", \"\" using "+k)+" ti \""+contig+"\"");
				++k;
				}
			makefileWriter.println("' | gnuplot");
			}
		
		
		if(!this.sample2stats.isEmpty())	
			{
			final String filename=toTsv("sample2gtype");
			PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
			pw.println("Sample\t"+Arrays.stream(GenotypeType.values()).map(T->T.name()).collect(Collectors.joining("\t")));
			for(final String sample: this.sample2stats.keySet())
				{
				pw.print(sample);
				for(final GenotypeType gtype: GenotypeType.values())
					{
					pw.print("\t"+this.sample2stats.get(sample).countTypes.count(gtype));
					}
				pw.println();
				}
			pw.flush();
			pw.close();
			
			
			final String png=toPng(filename);
			makefileWriter.println("ALL_TARGETS+=" + png);
			makefileWriter.println(png+":"+filename);

			makefileWriter.print("\techo 'set terminal png truecolor size ${SCREEN_WIDTH}, ${SCREEN_HEIGHT};"
					+ "set title \"Genotypes Types\";"
					+ "set xlabel \"Sample\";"
					+ "set key autotitle columnhead;"
					+ "set xtic rotate by 90 right;"
					+ "set ylabel \"Count\";"
					+ "set yrange [0:];"
					+ "set key invert reverse Left outside;"
					+ "set datafile separator \"\t\";"
					+ "set style fill solid border -1;"
					+ "set style data histograms;set style histogram rowstacked;"
					+ "set boxwidth 0.95;set output \"$@\";"
					+ "plot \"$<\" using 2:xtic(1)");
			int k=2;
			for(final GenotypeType gtype: GenotypeType.values())
				{
				makefileWriter.print((k==2?"":", \"\" using "+k)+" ti \""+gtype.name()+"\"");
				++k;
				}
			makefileWriter.println("' | gnuplot");
			}
			
			if(!this.countAffectedSamples.isEmpty())
				{
				final String filename=toTsv("affectedSamples");
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final RangeOfIntegers.Range k: VcfStats.this.affectedTranches.getRanges())
					{
					long n=this.countAffectedSamples.count(k); if(n==0L) continue;
					pw.println(k.toString()+"/"+this.sample2stats.size()+"\t"+n);
					}
				pw.flush();
				pw.close();
				
				final String png= toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
				makefileWriter.println("\techo '"
						+ "set ylabel \"Number of Variants /" + this.countVariants+"\";"
						+ "set xlabel \"Number of Affected Samples\";"
						+ "set xtic rotate by 90 right;"
						+ "set size  ratio 0.618;"
						+ "set title \"Number Variants/ num( Affected Samples ) \";"
						+ "set style fill solid border -1;"
						+ "set key  off;"
						+ "set datafile separator \"\t\";"
						+ "set auto x;"
						+ "set yrange [0:];"
						+ "set style histogram;"
						+ "set style data histogram;"
						+ "set terminal png truecolor;"
						+ "set output \"$@\";"
						+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
						+ "gnuplot");				
				}
				
			if(!this.countAltAlleles.isEmpty())
				{
				final String filename= toTsv("countAltAlleles");
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final RangeOfIntegers.Range k: VcfStats.this.altTranches.getRanges())
					{
					long n=this.countAltAlleles.count(k); if(n==0L) continue;
					pw.println(k.toString()+"\t"+n);
					}
				pw.flush();
				pw.close();
				
				final String png= toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
				makefileWriter.println("\techo '"
						+ "set ylabel \"Number of Variants /" + this.countVariants+"\";"
						+ "set xlabel \"Number of ALT alleles\";"
						+ "set size  ratio 0.618;"
						+ "set title \"Number of ALT alleles\";"
						+ "set style fill solid border -1;"
						+ "set key  off;"
						+ "set datafile separator \"\t\";"
						+ "set auto x;"
						+ "set yrange [0:];"
						+ "set style histogram;"
						+ "set style data histogram;"
						+ "set terminal png truecolor;"
						+ "set output \"$@\";"
						+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
						+ "gnuplot");				

				}		
			

			if(!this.countIndelSize.isEmpty())
				{
				final String filename= toTsv("countIndelSize");
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final RangeOfIntegers.Range k: VcfStats.this.indelTranches.getRanges())
					{
					long n=this.countIndelSize.count(k); if(n==0L) continue;
					pw.println(k.toString()+"\t"+n);
					}
				pw.flush();
				pw.close();
				
				final String png= toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
				makefileWriter.println("\techo '"
						+ "set ylabel \"Number of Variants /" + this.countVariants+"\";"
						+ "set yrange [0:];"
						+ "set xlabel \"Indel Size\";"
						+ "set xtic rotate by 90 right;"
						+ "set size  ratio 0.618;"
						+ "set title \"Indel Sizes\";"
						+ "set style fill solid border -1;"
						+ "set key  off;"
						+ "set datafile separator \"\t\";"
						+ "set auto x;"
						+ "set style histogram;"
						+ "set style data histogram;"
						+ "set terminal png truecolor;"
						+ "set output \"$@\";"
						+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
						+ "gnuplot");				

				}
			
			if(!this.countDepth.isEmpty())
				{
					{
					final String filename=toTsv("countDepth");
					PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
					for(final RangeOfIntegers.Range k: VcfStats.this.depthTranches.getRanges())
						{
						long n=this.countDepth.count(k); if(n==0L) continue;
						pw.println(k.toString()+"\t"+n);
						}
					pw.flush();
					pw.close();
					
					final String png= toPng(filename);
					makefileWriter.println("ALL_TARGETS+=" + png);
					makefileWriter.println(png+":"+filename);
					makefileWriter.println("\techo '"
							+ "set title \"Depth / Variant\";"
							+ "set ylabel \"Count / " + this.countVariants+" Variants\";"
							+ "set yrange [0:];"
							+ "set xlabel \"Depth\";"
							+ "set xtic rotate by 90 right;"
							+ "set size  ratio 0.618;"
							+ "set ylabel \"Depth /" + this.countVariants+" Variants\";"
							+ "set style fill solid border -1;"
							+ "set key  off;"
							+ "set datafile separator \"\t\";"
							+ "set auto x;"
							+ "set style histogram;"
							+ "set style data histogram;"
							+ "set terminal png truecolor;"
							+ "set output \"$@\";"
							+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
							+ "gnuplot");			
					}
				if(!this.sample2stats.isEmpty() && this.sample2stats.values().stream().filter(st->!st.countDepth.isEmpty()).count()>0) {

					final String filename=toTsv("countDepthBySample");
					final PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);

					pw.println("Sample\t"+ VcfStats.this.depthTranches.getRanges().stream().map(D->D.toString()).collect(Collectors.joining("\t")));
					for(final SampleStat st:this.sample2stats.values())
						{
						if(st.countDepth.isEmpty()) continue;//all sample are HOm_REF
						pw.print(st.sampleName);
						for(final RangeOfIntegers.Range range: VcfStats.this.depthTranches.getRanges())
							{
							pw.print("\t"+st.countDepth.count(range));
							}
						pw.println();
						}
					pw.flush();
					pw.close();
					
					
					final String png= toPng(filename);
					makefileWriter.println("ALL_TARGETS+=" + png);
					makefileWriter.println(png+":"+filename);
					makefileWriter.print("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT};"
							+ "set title \"Depth / Genotype\";"
							+ "set xlabel \"Sample\";"
							+ "set ylabel \"Depth Genotype\";"
							+ "set key autotitle columnhead;"
							+ "set xtic rotate by 90 right;"
							+ "set ylabel \"Count\";"
							+ "set yrange [0:];"
							+ "set key invert reverse Left outside;"
							+ "set datafile separator \"\t\";"
							+ "set style fill solid border -1;"
							+ "set style data histograms;set style histogram rowstacked;"
							+ "set boxwidth 0.95;set output \"$@\";"
							+ "plot \"$<\" using 2:xtic(1)"
							);
					int k=0;
					for(@SuppressWarnings("unused")final RangeOfIntegers.Range range: VcfStats.this.depthTranches.getRanges())
						{
						if(k>0) makefileWriter.print(",\"\" using "+(k+2));
						++k;
						}
		
					
					makefileWriter.println("'| gnuplot");
					
				}
				}
			if(!this.countDistances.isEmpty())
				{
					{
					final String filename=toTsv("countDistances");
					PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
					for(final RangeOfIntegers.Range k: VcfStats.this.distanceTranches.getRanges())
						{
						long n=this.countDistances.count(k); if(n==0L) continue;
						pw.println(k.toString()+"\t"+n);
						}
					pw.flush();
					pw.close();
					
					final String png= toPng(filename);
					makefileWriter.println("ALL_TARGETS+=" + png);
					makefileWriter.println(png+":"+filename);
					makefileWriter.println("\techo '"
							+ "set ylabel \"Count Variants\";"
							+ "set yrange [0:];"
							+ "set xlabel \"Distance\";"
							+ "set xtic rotate by 90 right;"
							+ "set size  ratio 0.618;"
							+ "set title \"Distance Between Variants\";"
							+ "set style fill solid border -1;"
							+ "set key  off;"
							+ "set datafile separator \"\t\";"
							+ "set auto x;"
							+ "set style histogram;"
							+ "set style data histogram;"
							+ "set terminal png truecolor;"
							+ "set output \"$@\";"
							+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
							+ "gnuplot");	
					}
				if(!this.sample2stats.isEmpty())
					{
					final String filename=toTsv("countDistancesBySample");
					final PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);

					pw.println("Sample\t"+ VcfStats.this.distanceTranches.getRanges().stream().map(D->D.toString()).collect(Collectors.joining("\t")));
					for(final SampleStat st:this.sample2stats.values())
						{
						if(st.countDistances.isEmpty()) continue;//all sample are HOm_REF
						pw.print(st.sampleName);
						for(final RangeOfIntegers.Range range: VcfStats.this.distanceTranches.getRanges())
							{
							pw.print("\t"+st.countDistances.count(range));
							}
						pw.println();
						}
					pw.flush();
					pw.close();
					
					
					final String png= toPng(filename);
					makefileWriter.println("ALL_TARGETS+=" + png);
					makefileWriter.println(png+":"+filename);
					makefileWriter.print("\techo 'set terminal png truecolor size ${SCREEN_WIDTH}, ${SCREEN_HEIGHT};"
							+ "set title \"Distance between variants\";"
							+ "set xlabel \"Sample\";"
							+ "set ylabel \"Distance between variants\";"
							+ "set key autotitle columnhead;"
							+ "set xtic rotate by 90 right;"
							+ "set ylabel \"Count\";"
							+ "set yrange [0:];"
							+ "set key invert reverse Left outside;"
							+ "set datafile separator \"\t\";"
							+ "set style fill solid border -1;"
							+ "set style data histograms;set style histogram rowstacked;"
							+ "set boxwidth 0.95;set output \"$@\";"
							+ "plot \"$<\" using 2:xtic(1)"
							);
					int k=0;
					for(@SuppressWarnings("unused")final RangeOfIntegers.Range range: VcfStats.this.distanceTranches.getRanges())
						{
						if(k>0) makefileWriter.print(",\"\" using "+(k+2));
						++k;
						}
		
					
					makefileWriter.println("'| gnuplot");
					}
				
				}
			
			if(this.mafPlotter!=null)
				{
				this.mafPlotter.close();
				
				final String png= toPng(this.mafPlotter.filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+this.mafPlotter.filename);
				makefileWriter.println("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT};"
						+ "set title \"MAF Cases/Controls\";"
						+ "set ylabel \"Controls\";"
						+ "set xlabel \"Cases\";"
						//+ "set style fill  transparent solid 0.35 noborder;"//http://stackoverflow.com/questions/34532568
						//+ "set style circle radius 0.02;" 
						+ "set nokey;"
						+ "set xrange [0:1];"
						+ "set yrange [0:1];"
						+ "set output \"$@\";"
						+ "plot \"$<\" u 1:2 with points;' | gnuplot");
				}
			
			if(!this.nucleicAcidChanges.isEmpty())
				{
				final String filename = toTsv("transvers");
				final PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				
				
				//pw.println("TYPE\tALL\tCDS");
				//pw.println("TRANSITION\t"+this.nucleicAcidChanges.count(NucleicAcidChange.transition)+"\t"+this.nucleicAcidChanges.count(NucleicAcidChange.transition_in_cds));
				//pw.println("TRANSVERSION\t"+this.nucleicAcidChanges.count(NucleicAcidChange.transversion)+"\t"+this.nucleicAcidChanges.count(NucleicAcidChange.transversion_in_cds));
				pw.println("TYPE\tTRANSITION\tTRANSVERSION");
				pw.println("ALL\t"+this.nucleicAcidChanges.count(NucleicAcidChange.transition)+"\t"+this.nucleicAcidChanges.count(NucleicAcidChange.transversion));
				pw.println("CDS\t"+this.nucleicAcidChanges.count(NucleicAcidChange.transition_in_cds)+"\t"+this.nucleicAcidChanges.count(NucleicAcidChange.transversion_in_cds));
				pw.flush();
				pw.close();
				
				final String png= toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
				makefileWriter.println("\techo 'set terminal png truecolor size ${SCREEN_WIDTH}, ${SCREEN_HEIGHT};"
						+ "set title \"Transition/Transversion\";"
						+ "set xlabel \"Where\";"
						+ "set ylabel \"Count Variants / "+this.countVariants+"\";"
						+ "set key autotitle columnhead;"
						+ "set xtic rotate by 90 right;"
						+ "set ylabel \"Count\";"
						+ "set yrange [0:];"
						+ "set key invert reverse Left outside;"
						+ "set datafile separator \"\t\";"
						+ "set style fill solid border -1;"
						+ "set style data histograms;set style histogram rowstacked;"
						+ "set boxwidth 0.95;set output \"$@\";"
						+ "plot \"$<\" using 2:xtic(1) , \"\" using 3 "
						+ "' | gnuplot");
				}
			
			if(!this.geneLocations.isEmpty())
				{
				final String filename = toTsv("geneLoc");
				final PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(GeneLocation loc:GeneLocation.values())
					{
					if(loc==GeneLocation.not_in_gene) continue;
					pw.println(loc.name()+"\t"+this.geneLocations.count(loc));
					}
				pw.flush();
				pw.close();
				
				
				final String png= toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
				makefileWriter.println("\techo '"
						+ "set ylabel \"Count Variants / "+this.countVariants+"\";"
						+ "set yrange [0:];"
						+ "set xlabel \"Category\";"
						+ "set xtic rotate by 90 right;"
						+ "set size  ratio 0.618;"
						+ "set title \"Position in the transcript\";"
						+ "set style fill solid border -1;"
						+ "set key  off;"
						+ "set datafile separator \"\t\";"
						+ "set auto x;"
						+ "set style histogram;"
						+ "set style data histogram;"
						+ "set terminal png truecolor;"
						+ "set output \"$@\";"
						+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
						+ "gnuplot");		

				}
			if(!this.consequences.isEmpty())
				{
				{
				final String filename=toTsv("predictions");
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final String k: this.consequences.keySet())
					{
				
					long n=this.consequences.count(k); if(n==0L) continue;
					pw.println(k.toString()+"\t"+n);
					}
				pw.flush();
				pw.close();
				
				final String png= toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
				makefileWriter.println("\techo '"
						+ "set ylabel \"Count Variants / "+this.countVariants+"\";"
						+ "set yrange [0:];"
						+ "set xlabel \"Consequences (warning: categories may overlap)\";"
						+ "set xtic rotate by 90 right;"
						+ "set size  ratio 0.618;"
						+ "set title \"Consequences/Variants\";"
						+ "set style fill solid border -1;"
						+ "set key  off;"
						+ "set datafile separator \"\t\";"
						+ "set auto x;"
						+ "set style histogram;"
						+ "set style data histogram;"
						+ "set terminal png truecolor;"
						+ "set output \"$@\";"
						+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
						+ "gnuplot");	
				}
				
				if(!this.sample2stats.isEmpty() && !this.sequenceOntologyTermsToObserve.isEmpty())
					{

					final String filename=toTsv("predictionsBySample");
					final PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);

					pw.println("Sample\t"+ this.sequenceOntologyTermsToObserve.stream().map(D->D.getLabel()).collect(Collectors.joining("\t")));
					for(final SampleStat st:this.sample2stats.values())
						{
						if(st.consequences.isEmpty()) continue;//all sample are HOm_REF
						pw.print(st.sampleName);
						for(final SequenceOntologyTree.Term term: this.sequenceOntologyTermsToObserve)
							{
							pw.print("\t"+st.consequences.count(term.getLabel()));
							}
						pw.println();
						}
					pw.flush();
					pw.close();
					
					
					final String png= toPng(filename);
					makefileWriter.println("ALL_TARGETS+=" + png);
					makefileWriter.println(png+":"+filename);
					makefileWriter.print("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT};"
							+ "set title \"Predictions/Samples\";"
							+ "set xlabel \"Sample\";"
							+ "set ylabel \"Prediction\";"
							+ "set key autotitle columnhead;"
							+ "set xtic rotate by 90 right;"
							+ "set ylabel \"Count\";"
							+ "set yrange [0:];"
							+ "set key invert reverse Left outside;"
							+ "set datafile separator \"\t\";"
							+ "set style fill solid border -1;"
							+ "set style data histograms;set style histogram rowstacked;"
							+ "set boxwidth 0.95;set output \"$@\";"
							+ "plot \"$<\" using 2:xtic(1)"
							);
					int k=0;
					for(@SuppressWarnings("unused") final SequenceOntologyTree.Term term: this.sequenceOntologyTermsToObserve)
						{
						if(k>0) makefileWriter.print(",\"\" using "+(k+2));
						++k;
						}
		
					
					makefileWriter.println("'| gnuplot");
					}
				
				//mendelian
				if(!this.sample2stats.isEmpty() && VcfStats.this.pedigree.hasTrios())
					{
					final String filename=toTsv("mendel");
					final PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);

					pw.println("Sample\t"+ this.sequenceOntologyTermsToObserve.stream().map(D->D.getLabel()).collect(Collectors.joining("\t")));
					for(final SampleStat st:this.sample2stats.values())
						{
						pw.print(st.sampleName);
						for(final SequenceOntologyTree.Term term: this.sequenceOntologyTermsToObserve)
							{
							pw.print("\t"+st.countMendelianViolations.count(term.getLabel()));
							}
						pw.println();
						}
					pw.flush();
					pw.close();
					
					
					final String png= toPng(filename);
					makefileWriter.println("ALL_TARGETS+=" + png);
					makefileWriter.println(png+":"+filename);
					makefileWriter.print("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT};"
							+ "set title \"Mendelian violations\";"
							+ "set xlabel \"Sample\";"
							+ "set ylabel \"Prediction\";"
							+ "set key autotitle columnhead;"
							+ "set xtic rotate by 90 right;"
							+ "set ylabel \"Count\";"
							+ "set yrange [0:];"
							+ "set key invert reverse Left outside;"
							+ "set datafile separator \"\t\";"
							+ "set style fill solid border -1;"
							+ "set style data histograms;"
							+ "set style histogram rowstacked;"
							+ "set boxwidth 0.95;set output \"$@\";"
							+ "plot \"$<\" using 2:xtic(1)"
							);
					int k=0;
					for(@SuppressWarnings("unused")final SequenceOntologyTree.Term term: this.sequenceOntologyTermsToObserve)
						{
						if(k>0) makefileWriter.print(",\"\" using "+(k+2));
						++k;
						}
		
					makefileWriter.println("'| gnuplot");
					}

				}
			
			if(!VcfStats.this.disableGenotypeConcordance && !this.genotypeConcordance.isEmpty())
				{
				final String filename = toTsv("gtConcordance");
				final PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				
				//http://gnuplot.sourceforge.net/demo/heatmaps.html
				for(int x=0;x<  VcfStats.this.sampleNamesInOrder.size();++x)
					{
					pw.print(",");
					pw.print( VcfStats.this.sampleNamesInOrder.get(x));
					}
				pw.println();
				for(int y=0;y<  VcfStats.this.sampleNamesInOrder.size();++y)
					{
					pw.print( VcfStats.this.sampleNamesInOrder.get(y));
					for(int x=0;x<  VcfStats.this.sampleNamesInOrder.size();++x)
						{
						pw.print(",");
						pw.print(this.genotypeConcordance.count(new SamplePair(x,y)));
						}
					pw.println();
					}
				
				pw.flush();
				pw.close();
				
				
				final String png= toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
				makefileWriter.print("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT};"
						+ "set title \"Genotype Concordance\";"
						+ "unset key; set view map;"
						+ "set output \"$@\";"
						+ "set datafile separator \",\";"
						+ "plot \"$<\" matrix rowheaders columnheaders using 1:2:3 with image pixels;"
						);
				makefileWriter.println("'| gnuplot");
				}
			if(!this.countStructuralVariations.isEmpty())
				{
				final String filename= toTsv("countSV");
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final StructuralVariantType k: StructuralVariantType.values())
					{
					final long n=this.countStructuralVariations.count(k);
					pw.println(k.name()+"\t"+n);
					}
				pw.flush();
				pw.close();
				
				final String png= toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
				makefileWriter.println("\techo '"
						+ "set ylabel \"Number of Variants /" + this.countVariants+"\";"
						+ "set yrange [0:];"
						+ "set xlabel \"Structural Variation type\";"
						+ "set xtic rotate by 90 right;"
						+ "set size  ratio 0.618;"
						+ "set title \"Structural Variations\";"
						+ "set style fill solid border -1;"
						+ "set key  off;"
						+ "set datafile separator \"\t\";"
						+ "set auto x;"
						+ "set style histogram;"
						+ "set style data histogram;"
						+ "set terminal png truecolor;"
						+ "set output \"$@\";"
						+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
						+ "gnuplot");				
	
				}
			
			if(!this.sample2stats.values().stream().anyMatch(P->P.countStructuralVariations.isEmpty()))
				{
				final String filename=toTsv("sample2sv");
				
				
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				pw.println("Sample\t"+Arrays.stream(StructuralVariantType.values()).map(V->V.name()).collect(Collectors.joining("\t")));
				for(final String sample: this.sample2stats.keySet())
					{
					pw.print(sample);
					for(final StructuralVariantType svt: StructuralVariantType.values())
						{
						pw.print("\t"+this.sample2stats.get(sample).countStructuralVariations.count(svt));
						}
					pw.println();
					}
				pw.flush();
				pw.close();
				
				final String png=toPng(filename);
				makefileWriter.println("ALL_TARGETS+=" + png);
				makefileWriter.println(png+":"+filename);
	
				makefileWriter.print("\techo 'set terminal png truecolor size ${SCREEN_WIDTH},${SCREEN_HEIGHT} ;"
						+ "set title \"Sample/Structural Variation\";"
						+ "set xlabel \"Sample\";"
						+ "set key autotitle columnhead;"
						+ "set xtic rotate by 90 right;"
						+ "set ylabel \"Count\";"
						+ "set yrange [0:];"
						+ "set key invert reverse Left outside;"
						+ "set datafile separator \"\t\";"
						+ "set style fill solid border -1;"
						+ "set style data histograms;set style histogram rowstacked;"
						+ "set boxwidth 0.95;set output \"$@\";"
						+ "plot \"$<\" using 2:xtic(1)");
				int k=2;
				for(final StructuralVariantType svt: StructuralVariantType.values())
					{
					makefileWriter.print((k==2?"":", \"\" using "+k)+" ti \""+svt.name()+"\"");
					++k;
					}
				makefileWriter.println("' | gnuplot");
				}
			
			
			for(final SampleStat st:this.sample2stats.values())
				{
				st.finish(makefileWriter);
				}
			}
		
		}
	
	public VcfStats()
		{
		//this.selectExpressions.add("vc azd");
		}
	
	public List<KnownGene> getOverlappingKnownGenes(final VariantContext ctx)
		{
		if(this.knownGeneTreeMap==null) return Collections.emptyList();
		final List<KnownGene> L = new ArrayList<>();
		for(final List<KnownGene> lkg:VcfStats.this.knownGeneTreeMap.getOverlapping(ctx))
			{
			L.addAll(lkg);
			}
		return L;
		}
	
	// https://en.wikipedia.org/wiki/File:Transitions-transversions-v3.png
	private static boolean isTransversion(final Character a1, final Character a2)
		{
		if(a1==null || a2==null) return false;
		if(a1=='A' &&  a2=='C') return true;
		if(a1=='C' &&  a2=='A') return true;
		if(a1=='G' &&  a2=='T') return true;
		if(a1=='T' &&  a2=='G') return true;
		return false;
		}

	private static boolean isTransition(final Character a1, final Character a2)
		{
		if(a1==null || a2==null) return false;
		if(a1=='A' &&  a2=='G') return true;
		if(a1=='G' &&  a2=='A') return true;
		if(a1=='C' &&  a2=='T') return true;
		if(a1=='T' &&  a2=='C') return true;
		return false;
		}

	private static Character asSimpleATGC(final Allele al)
		{
		if(al==null || al.isSymbolic() || al.equals(Allele.SPAN_DEL)) return null;
		final String s=al.getBaseString().toUpperCase();
		if(s==null || s.equals(".") || s.length()!=1 ) return null;
		switch(s.charAt(0))
			{
			case 'A': case 'T': case 'G': case 'C': return s.charAt(0);
			default: return null;
			}
		}

	
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.binSize<=0) {
			LOG.error("binSize < 0");
			return -1;
		}
		
		VariantContextWriter teeOut=null;
		VCFIterator iter = null;
		final Map<String,VariantStats> category2stats = new HashMap<>();
		
		PrintWriter makefileWriter =null;
		try {
			
			
			this.archiveFactory = ArchiveFactory.open(this.outputFile);
			if(this.tee) teeOut = super.openVariantContextWriter(null);
			
			iter= super.openVCFIterator(oneFileOrNull(args));
			
			
			
			final VCFHeader header=iter.getHeader();
			this.sampleNamesInOrder = Collections.unmodifiableList(header.getSampleNamesInOrder());
			
			final SAMSequenceDictionary dict=header.getSequenceDictionary();
			if(dict!=null && !dict.isEmpty()) {
				this.the_dictionary = dict;
				}
			
			if(this.kgFile!=null)
				{
				LOG.info("load "+kgFile);
				this.knownGeneTreeMap=KnownGene.loadUriAsIntervalTreeMap(this.kgFile,KG->(dict==null || dict.getSequence(KG.getContig())!=null));
				}
			else
				{
				this.knownGeneTreeMap=null;
				}
			if(this.pedigreeFile!=null)
				{
				this.pedigree = Pedigree.newParser().parse(this.pedigreeFile);
				}
			else
				{
				Pedigree tmpPed=null;
				try 
					{
					tmpPed =  Pedigree.newParser().parse(header);
					}
				catch(Exception err) {
					tmpPed = Pedigree.createEmptyPedigree();
					}
				this.pedigree = tmpPed;
				}
			makefileWriter = this.archiveFactory.openWriter(this.prefix+"Makefile");
			makefileWriter.println(".PHONY: all all_targets ");
			makefileWriter.println("SCREEN_WIDTH?=2600");
			makefileWriter.println("SCREEN_HEIGHT?=1000");
			makefileWriter.println("ALL_TARGETS=");
			makefileWriter.println("all: all_targets");


			
			
			if(teeOut!=null) teeOut.writeHeader(header);
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header).logger(LOG);
			while(iter.hasNext())
				{
				final VariantContext ctx=progress.watch(iter.next());
				if(teeOut!=null) teeOut.add(ctx);

				
				for(final String category: this.variantToCategoryKeys.apply(ctx))
					{
					VariantStats vcstat = category2stats.get(category);
					if(vcstat==null) {
						vcstat = new VariantStats(category,header);
						category2stats.put(category, vcstat);
						}
					vcstat.visit(ctx);
					}
				
				
				}
			for(final String category: category2stats.keySet())
				{	
				final VariantStats vcstats = category2stats.get(category);
				vcstats.finish(makefileWriter);
				}

			progress.finish();
			makefileWriter.println("all_targets : ${ALL_TARGETS}");
			makefileWriter.flush();makefileWriter.close();makefileWriter=null;
			
			iter.close();iter=null;
			this.archiveFactory.close();archiveFactory=null;
			if(teeOut!=null) teeOut.close(); teeOut=null;
			return 0;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
		} finally
			{
			knownGeneTreeMap=null;
			CloserUtil.close(archiveFactory);
			CloserUtil.close(teeOut);
			CloserUtil.close(iter);
			CloserUtil.close(makefileWriter);
			}
		
		}
	
			
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new VcfStats().instanceMainWithExit(args);
		}
	}
