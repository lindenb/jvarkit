/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.math.RangeOfIntegers;
import com.github.lindenb.jvarkit.tools.burden.MafCalculator;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.semontology.Term;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
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
$ java -jar dist/vcfstats.jar ~karaka/BURDEN_JVARKIT/MVP/20170227.pct001gBED.Q.mvp_frex.vcf.gz -o tmp
$ ls tmp/
ALL.sample2gtype.tsv  ALL.variant2type.tsv  Makefile

$ head tmp/ALL.sample2gtype.tsv  tmp/ALL.variant2type.tsv

==> tmp/ALL.sample2gtype.tsv <==
Type	NO_CALL	HOM_REF	HET	HOM_VAR	UNAVAILABLE	MIXED
X0G73J	2538	3440	218	132	0	0
Y00G3I	2543	3462	193	130	0	0
Z03K	2547	3417	252	112	0	0
A0G3N	2547	3424	209	148	0	0
B980	1909	4068	202	149	0	0
C003P	2559	3417	216	136	0	0
D0741	2557	3428	204	139	0	0
E073O	2566	3433	212	117	0	0
F00G7	2560	3444	191	133	0	0
(...)

==> tmp/ALL.variant2type.tsv <==
Type	Count
MIXED	7
SNP	6125
INDEL	196

```


END_DOC
 */
@Program(name="vcfstats",
	description="Produce VCF statitics",
	keywords={"vcf","stats","burden","gnuplot"},
	terms=Term.ID_0000018
	)
public class VcfStats extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfStats.class).make();

	@Parameter(names={"-o","--output"},description="output Directory or zip file",required=true)
	private File outputFile = null;
	
	@Parameter(names={"-kg","--knownGenes"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private File kgFile = null;
	
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
	private RangeOfIntegers affectedTranches = new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,20,50,100);
	@Parameter(names={"--trancheDP"},description="tranches for the DEPTH. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers depthTranches =new RangeOfIntegers(0,10,20,30,50,100,200);
	@Parameter(names={"--trancheIndelSize"},description="tranches for the Indel size "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers indelTranches =new RangeOfIntegers(0,1,2,3,4,5,6,8,9,10,15,20);
	@Parameter(names={"--trancheAlts"},description="tranches for the number of ALTs. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers altTranches =new RangeOfIntegers(0,1,2,3,4,5,6,8,9,10);
	@Parameter(names={"--trancheDistance"},description="tranches for the distance between the variants. "+RangeOfIntegers.OPT_DESC,converter=RangeOfIntegers.StringConverter.class)
	private RangeOfIntegers distanceTranches =new RangeOfIntegers(0,1,2,3,4,5,6,7,8,9,10,20,100,200,300,400,500,1000);

	
	private ArchiveFactory archiveFactory=null;
	
	
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
	
	private final Function<VariantContext, Set<String>> variantToMafKeys = VC ->{
		if( mafKeyInINFO==null || mafKeyInINFO.trim().isEmpty()) {
			return Collections.singleton("ALL");
			}
		else
			{
			return VC.getAttributeAsList(mafKeyInINFO).
				stream().
				filter(O->!(O==null || ".".equals(O))).
				map(O->String.valueOf(O)).
				collect(Collectors.toSet());
			}
		};

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

		
	private abstract class AbstractStat
		{
		final Counter<RangeOfIntegers.Range> countDepth = new Counter<>();

		public abstract void visit(final VariantContext ctx);
		}
	
	private class SampleStat extends AbstractStat
		{
		final Counter<GenotypeType> countTypes = new Counter<>();

		final String sampleName;
		SampleStat(final String sampleName) {
			this.sampleName = sampleName;
			}
		@Override
		public void visit(final VariantContext ctx) {
			final Genotype genotype = ctx.getGenotype(this.sampleName);
			if(genotype==null) return;
			this.countTypes.incr(genotype.getType());
			
			
			if(genotype.hasDP())
				{
				final int dp = genotype.getDP();
				if(dp>=0)
					{
					this.countDepth.incr(VcfStats.this.depthTranches.getRange(dp));
					}
				}
			}
		}
	
	private class VariantStats extends AbstractStat
		{
		private final String key;
		private PlotMaf mafPlotter= null;
		private final Set<String> affectedSamples;
		private final Set<String> unaffectedSamples;
		private final Map<String,SampleStat> sample2stats = new TreeMap<>();
		final Counter<VariantContext.Type> countTypes = new Counter<>();
		final Counter<RangeOfIntegers.Range> countAffectedSamples = new Counter<>();
		final Counter<RangeOfIntegers.Range> countAltAlleles = new Counter<>();
		final Counter<RangeOfIntegers.Range> countIndelSize = new Counter<>();
		final Counter<RangeOfIntegers.Range> countDistances = new Counter<>();
		int count_transition=0;
		int count_transversion=0;
		private ContigPosRef prevCtx=null;
		
		VariantStats(final String key,final VCFHeader header) {
			this.key = key;
			for(final String sn: header.getSampleNamesInOrder())
				{
				this.sample2stats.put(sn,new SampleStat(sn));
				}
			
			this.affectedSamples = 
					VcfStats.this.pedigree.getPersons().stream().
						filter(P->P.isAffected()).
						map(P->P.getId()).
						filter(S->header.getSampleNamesInOrder().contains(S)).
						collect(Collectors.toSet())
						;
			
			this.unaffectedSamples = 
					VcfStats.this.pedigree.getPersons().stream().
						filter(P->P.isUnaffected()).
						map(P->P.getId()).
						filter(S->header.getSampleNamesInOrder().contains(S)).
						collect(Collectors.toSet())
						;

			}
		
		public void visit(final VariantContext ctx) {
			
			this.countTypes.incr(ctx.getType());
			for(final SampleStat st: this.sample2stats.values()) st.visit(ctx);
			
			//distance
			final ContigPosRef contigPosRef = new ContigPosRef(ctx);
			if(prevCtx!=null && prevCtx.getContig().equals(contigPosRef.getContig()) && prevCtx.getStart() <= contigPosRef.getStart())
				{
				final int distance = contigPosRef.getStart() - this.prevCtx.getStart();
				this.countDistances.incr(VcfStats.this.distanceTranches.getRange(distance));
				}
			prevCtx=contigPosRef;
			
			// MAF
			
			final List<Allele> alternates = ctx.getAlternateAlleles();

			if(alternates.size()==1 )
				{
				if(isTransition(asSimpleATGC(ctx.getReference()),asSimpleATGC(alternates.get(0)))) {
					this.count_transition++;
					}	
				if(isTransversion(asSimpleATGC(ctx.getReference()),asSimpleATGC(alternates.get(0)))) {
					this.count_transversion++;
					}	
				}
			
			/* can we compute the MAF ? we need (non)affected samples*/
			if(!(alternates.isEmpty() || this.unaffectedSamples.isEmpty() || this.affectedSamples.isEmpty()))
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
							filter(G->!(G.isNoCall() || G.isHomRef() || G.isFiltered() )).
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
			if(alternates.size()>1)
				{
				this.countAltAlleles.incr(VcfStats.this.altTranches.getRange(alternates.size()));
				}
			

			
			}
		
		void finish(PrintWriter makefileWriter) throws IOException
		{

			{
			final String filename=prefix+this.key+".variant2type.tsv";

			PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
			pw.println("Type\tCount");

			for(final VariantContext.Type type: this.countTypes.keySet())
				{
				pw.println(type.name()+"\t"+this.countTypes.count(type));
				}
			pw.flush();
			pw.close();
			final String png= "$(patsubst %.tsv,%.png,"+filename+")";
			makefileWriter.println("ALL_TARGETS+=" + png);
			makefileWriter.println(png+":"+filename);
			makefileWriter.println("\techo 'set key autotitle columnheader;set ylabel \"Count\";set size  ratio 0.618;set title \"Variant Types\";set style fill solid border -1;set key  off;set datafile separator \"\t\";set auto x;set style histogram;set style data histogram;set terminal png truecolor  ;set output \"$@\";plot \"$<\" using 2:xtic(1) ti \"Variant Type\";' | gnuplot");
			}

			{
			final String filename=prefix+this.key+".sample2gtype.tsv";
			PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
			pw.println("Type\t"+Arrays.stream(GenotypeType.values()).map(T->T.name()).collect(Collectors.joining("\t")));
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
			
			
			final String png= "$(patsubst %.tsv,%.png,"+filename+")";
			makefileWriter.println("ALL_TARGETS+=" + png);
			makefileWriter.println(png+":"+filename);

			makefileWriter.print("\techo 'set terminal png truecolor size 2600, 1000;"
					+ "set title \"Genotypes Types\";set xlabel \"Sample\";"
					+ "set xtic rotate by 90;set ylabel \"Count\";"
					+ "set key invert reverse Left outside;set datafile separator \"\t\";"
					+ "set style fill solid border -1;set style data histograms;set style histogram rowstacked;"
					+ "set boxwidth 0.95;set output \"$@\";plot \"$<\" using 2:xtic(1)");
			int k=0;
			for(final GenotypeType gtype: GenotypeType.values())
				{
				makefileWriter.print((k==0?"":", \"\" using "+k)+" ti \""+gtype.name()+"\"");
				++k;
				}
			makefileWriter.println("' | gnuplot");
			}
			
			if(!this.countAffectedSamples.isEmpty())
				{
				final String filename=prefix+this.key+".affectedSamples.tsv";
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final RangeOfIntegers.Range k: VcfStats.this.affectedTranches.getRanges())
					{
					pw.println(k.toString()+"\t"+this.countAffectedSamples.count(k));
					}
				pw.flush();
				pw.close();
				}
				
			if(!this.countAltAlleles.isEmpty())
				{
				final String filename=prefix+this.key+".countAltAlleles.tsv";
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final RangeOfIntegers.Range k: VcfStats.this.altTranches.getRanges())
					{
					pw.println(k.toString()+"\t"+this.countAltAlleles.count(k));
					}
				pw.flush();
				pw.close();
				}	
			if(!this.countIndelSize.isEmpty())
				{
				final String filename=prefix+this.key+".countIndelSize.tsv";
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final RangeOfIntegers.Range k: VcfStats.this.indelTranches.getRanges())
					{
					pw.println(k.toString()+"\t"+this.countIndelSize.count(k));
					}
				pw.flush();
				pw.close();
				}
			if(!this.countDepth.isEmpty())
				{
				final String filename=prefix+this.key+".countDepth.tsv";
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final RangeOfIntegers.Range k: VcfStats.this.depthTranches.getRanges())
					{
					pw.println(k.toString()+"\t"+this.countDepth.count(k));
					}
				pw.flush();
				pw.close();
				}
			if(!this.countDistances.isEmpty())
				{
				final String filename=prefix+this.key+".countDistances.tsv";
				PrintWriter pw = VcfStats.this.archiveFactory.openWriter(filename);
				for(final RangeOfIntegers.Range k: VcfStats.this.distanceTranches.getRanges())
					{
					pw.println(k.toString()+"\t"+this.countDistances.count(k));
					}
				pw.flush();
				pw.close();
				}
			
			
		}
		
		}
	
	public VcfStats()
		{
		//this.selectExpressions.add("vc azd");
		}
	
	// https://en.wikipedia.org/wiki/File:Transitions-transversions-v3.png
	private static boolean isTransversion(Character a1, Character a2)
		{
		if(a1==null || a2==null) return false;
		if(a1=='A' &&  a2=='C') return true;
		if(a1=='C' &&  a2=='A') return true;
		if(a1=='G' &&  a2=='T') return true;
		if(a1=='T' &&  a2=='G') return true;
		return false;
		}

	private static boolean isTransition(Character a1, Character a2)
		{
		if(a1==null || a2==null) return false;
		if(a1=='A' &&  a2=='G') return true;
		if(a1=='G' &&  a2=='A') return true;
		if(a1=='C' &&  a2=='T') return true;
		if(a1=='T' &&  a2=='C') return true;
		return false;
		}

	
	
	private static Character asSimpleATGC(Allele al)
		{
		if(al==null) return null;
		String s=al.getBaseString().toUpperCase();
		if(s==null || s.equals(".") || s.length()!=1 ) return null;
		switch(s.charAt(0))
			{
			case 'A': case 'T': case 'G': case 'C': return s.charAt(0);
			default: return null;
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		VariantContextWriter teeOut=null;
		try {
			
			
			this.archiveFactory = ArchiveFactory.open(this.outputFile);
			if(this.tee) teeOut = super.openVariantContextWriter(null);
			scan(oneFileOrNull(args),teeOut);
			archiveFactory.close();
			archiveFactory=null;
			if(teeOut!=null) teeOut.close();
			teeOut=null;
			return 0;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
		} finally
			{
			CloserUtil.close(archiveFactory);
			CloserUtil.close(teeOut);
			}
		
		}
	
	private void scan(final String vcfInputStream,VariantContextWriter teeOut) {
		VcfIterator iter = null;
		final Map<String,VariantStats> category2stats = new HashMap<>();
		
		PrintWriter makefileWriter =null;
		try
			{			
			iter= super.openVcfIterator(vcfInputStream);
			final VCFHeader header=iter.getHeader();
			
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
			makefileWriter.println("ALL_TARGETS=");
			makefileWriter.println("all: all_targets");

			
			final Set<String> intersectSamples = this.pedigree.getPersons().stream().
						map(P->P.getId()).
						filter(S->header.getSampleNamesInOrder().contains(S)).
						collect(Collectors.toSet())
						;
			


				/*
				generic_maf_gnuplot_filename = prefix+"maf.gnuplot~";
				PrintWriter pw= archiveFactory.openWriter(generic_maf_gnuplot_filename);
				pw.println("set title \"__TITLE__\" ; set ylabel \"Control\"; set xlabel \"Case\";");
				//http://stackoverflow.com/questions/34532568
				pw.println("set style fill  transparent solid 0.35 noborder");
				pw.println("set style circle radius 0.02");
				pw.println("set nokey ; set terminal  png truecolor ; ");
				pw.println("set xrange [0:1];");
				pw.println("set yrange [0:1];");
				pw.println("set output '__OUTPUT__';");
				pw.println("plot '__INPUT__' u 1:2 with circles lc rgb \"blue\"");
				pw.flush();
				pw.close();
				*/
				
			
			
			if(teeOut!=null) teeOut.writeHeader(header);
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
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
				if(header.hasGenotypingData())
					{
					final Set<String> variantContextKeys =  this.variantToMafKeys.apply(ctx);
					
				
	
					
					
					
					
					}
				
				
				}
			for(final String category: category2stats.keySet())
				{	
				final VariantStats vcstats = category2stats.get(category);
				vcstats.finish(makefileWriter);
				}

			
			
			progress.finish();
			iter.close();
			makefileWriter.println("all_targets : ${ALL_TARGETS}");
			makefileWriter.flush();makefileWriter.close();makefileWriter=null;
			iter=null;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			throw new RuntimeException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(vcfInputStream);
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
