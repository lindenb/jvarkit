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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.vcf.VCFHeader;

/**

BEGIN_DOC


### Synopsis

```
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) stdin 
$ java -jar dist/vcfcomparecallers.jar file1.vcf(.gz) file2.vcf(.gz) 

```


both vcf **must** share the same sequence dictionary and must be sorted

#### History

* 20170704 : rewritten from scratch

### Example


```
$ java -jar dist/vcfcomparecallers.jar -o tmp Proj1.samtools.vcf.gz  Proj1.varscan.vcf.gz
$ (cd tmp && make)
```


END_DOC
*/

@Program(name="vcfcomparecallers",
	description="Compare two VCFs and print common/exclusive information for each sample/genotype",
	keywords={"vcf","compare","genotype"}
	)
public class VcfCompareCallers
	extends Launcher
	{

	private static final Logger LOG = Logger.build(VcfCompareCallers.class).make();


	@Parameter(names={"-o","--output"},description="Directory or zip file to save results to be plotted with gnuplot",required=true)
	private File archiveFile = null;
	@Parameter(names={"-p","--prefix"},description="Archive prefix (for option -d)")
	private String archivePrefix ="";

	//@Parameter(names={"-n","--num"},description="number of variants to dump in the example file")
	//private int numberOfExampleVariants = 10 ;

	//@Parameter(names={"-e","--examplefile"},description="Write a few Variants in this XML file. Optional")
	//private File exampleFile = null;

	@Parameter(names={"-B","--bed"},description="Limit to variants in that BED region")
	private File captureFile = null;

	@Parameter(names={"-c","--noCallIsHomRef"},description="No Call is HomRef (created when comparing merged vcf with GATK: there is no homref, everything is nocall)")
	private boolean noCallIsHomRef = false;

	
	@Parameter(names={"-vcf1","--vcf1"},description="short descriptive name for VCF1")
	private String vcf1Name = "VCF1";
	@Parameter(names={"-vcf2","--vcf2"},description="short descriptive name for VCF2")
	private String vcf2Name = "VCF2";
	@Parameter(names={"--collapseGenotypeType"},description="collapse Genotype Type. Just show Same Genotype or Discordant, don't print the type of genotype.")
	private boolean collapseGenotypeType = false;
	@Parameter(names={"--jexl1"},description="An optional list of GATK-like JEXL expressions to filter the variants from VCF File 1")
	private List<String> jexlExprStrings1 = new ArrayList<>();
	@Parameter(names={"--jexl2"},description="An optional list of GATK-like JEXL expressions to filter the variants from VCF File 2")
	private List<String> jexlExprStrings2 = new ArrayList<>();
	

	// https://stackoverflow.com/questions/13655048/
	private static final String escapeUnderscore(final String s)
		{
		return s.replace('_', '-');
		}
	
	private class SampleCategory
		{
		final String variantCatName;
		final Counter<String> counter = new Counter<>();
		SampleCategory(final String variantCatName)
			{
			this.variantCatName=variantCatName;
			}
		void visit(final Genotype g1,final Genotype g2)
			{
			if(g1==null && g2==null)
				{
				counter.incr("Both Missing");
				}
			else if(g1==null && g2!=null)
				{
				counter.incr("Unique in "+vcf2Name);
				}
			else if(g1!=null && g2==null)
				{
				counter.incr("Unique in "+vcf1Name);
				}
			else
				{
				if(g1.sameGenotype(g2))
					{
					counter.incr(collapseGenotypeType?"Same":" "+g1.getType().name());
					}
				else
					{
					counter.incr(collapseGenotypeType?"Discordant":vcf1Name+" "+g1.getType().name()+ " -> "+g2.getType().name()+" "+vcf2Name);
					}
				}
			}
		}
	
	private class SampleInfo
		{
		final String sampleName;
		final Map<String,SampleCategory> sampleCat = new HashMap<>();
		SampleInfo(final String sampleName)
			{
			this.sampleName=sampleName;
			}
		
		/* return a non null variant from a list of two, check they're consistent */
		private VariantContext theOne(final VariantContext ctx0,final VariantContext ctx1)
			{
			if(ctx0==null && ctx1==null) throw new IllegalStateException();
			if(ctx0==null && ctx1!=null) return ctx1;
			if(ctx1==null && ctx0!=null) return ctx0;
			if(VcfCompareCallers.this.compareChromPosRef.compare(ctx0, ctx1)!=0) {
				throw new IllegalStateException();
				}
			return ctx0;
			}
		
		private Genotype makeHomRef(final VariantContext one,int ploidy)
			{
			final List<Allele> L = new ArrayList<>(ploidy);
			while(L.size()<ploidy) L.add(one.getReference());
			return new GenotypeBuilder(this.sampleName).alleles(L).make();
			}
		
		void visit(final VariantContext ctx0,final VariantContext ctx1)
			{
			final Set<String> categories = new HashSet<>();
			categories.add("ALL");
			final VariantContext theOne = theOne(ctx0,ctx1);
			categories.add(theOne.getType().name());
			
			
			Genotype g0=(ctx0==null?null:ctx0.getGenotype(this.sampleName));
			Genotype g1=(ctx1==null?null:ctx1.getGenotype(this.sampleName));
			
			if(g0!=null &&  VcfCompareCallers.this.noCallIsHomRef && !g0.isCalled())
				{
				g0 =  makeHomRef( this.theOne(ctx0, ctx1),g0.getPloidy());
				}
			if(g1!=null &&  VcfCompareCallers.this.noCallIsHomRef && !g1.isCalled())
				{
				g1 = makeHomRef( this.theOne(ctx0, ctx1),g1.getPloidy());
				}
			for(final String cat:categories)
				{
				SampleCategory sc = this.sampleCat.get(cat);
				if(sc==null) {
					sc = new SampleCategory(cat);
					 this.sampleCat.put(cat,sc);
					}
				sc.visit(g0,g1);
				}
			}
		}
	
	

	
	public VcfCompareCallers()
		{
		}
	
	private SAMSequenceDictionary global_dictionary=null;
	
	private final Function<String,Integer> contig2tid=C->{
		final int tid = global_dictionary.getSequenceIndex(C);
		if(tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(C, global_dictionary);
		return tid;
		};
	
	private final Comparator<String> compareContigs = (C1,C2)->{
		if(C1.equals(C2)) return 0;
		return contig2tid.apply(C1) - contig2tid.apply(C2);
		};
	
	private final Comparator<VariantContext> compareChromPos = (V1,V2)->{
		int i = compareContigs.compare(V1.getContig(),V2.getContig());
		if( i!=0 ) return i;
		return V1.getStart() - V2.getStart();
		};	
	private final Comparator<VariantContext> compareChromPosRef = (V1,V2)->{
		int i = compareChromPos.compare(V1,V2);
		if( i!=0 ) return i;
		return V1.getReference().compareTo(V2.getReference());
		};	
		

	/** container uri+vcfIterator */
	private class PeekVCF implements Closeable
		{
		final String uri;
		final VCFIterator iter;
		final VCFHeader header;
		final List<VariantContext> buffer = new ArrayList<>();
		final SAMSequenceDictionary dict;
		private int debug_multiple=0;
		int count=0;
		private final List<VariantContextUtils.JexlVCMatchExp> jexlVCMatchExps = new ArrayList<>();
		PeekVCF(final VCFIterator iterator,final String uri) throws IOException {
			this.uri = uri;
			this.iter = iterator;
			this.header = this.iter.getHeader();
			this.dict = this.header.getSequenceDictionary();
			if(this.dict==null || this.dict.isEmpty())
				{
				throw new JvarkitException.DictionaryMissing(uri);
				}
			}
		
		private List<VariantContext> priv_peek()
			{
			if(!this.buffer.isEmpty()) return this.buffer;
			while(this.iter.hasNext())
				{
				final VariantContext ctx= this.iter.peek();
				if(!ctx.isVariant())
					{
					this.iter.next();//consume
					continue;
					}
				
				if(!this.jexlVCMatchExps.isEmpty())
					{
					if(VariantContextUtils.match(ctx,this.jexlVCMatchExps).
						values().
						stream().
						anyMatch(B->B.booleanValue()==false)
						)
						{
						this.iter.next();//consume
						continue;
						}
					}
				
				if(this.buffer.isEmpty())
					{
					this.buffer.add(this.iter.next());
					this.count++;
					}
				else
					{
					// compare with first item in buffer
					final int i = VcfCompareCallers.this.compareChromPos.compare(
							ctx,
							this.buffer.get(0) 
							);
					if( i< 0) {
						throw new JvarkitException.UserError("Variant are not sorted! got: "+ctx+" after "+buffer.get(0));
						}
					else if(i > 0)
						{
						break;
						}
					else //i == 0 or growing
						{
						buffer.add(this.iter.next());
						this.count++;
						}
					}
				}
			Collections.sort(this.buffer, VcfCompareCallers.this.compareChromPosRef);
			return buffer;
			}
		VariantContext peek()
			{
			final List<VariantContext> L = priv_peek();
			if(L.isEmpty()) return null;
			if(L.size()==1) return L.get(0);
			
			final List<VariantContext> L2 = L.stream().
					filter(V->V==L.get(0) /* compare ptr */|| VcfCompareCallers.this.compareChromPosRef.compare(L.get(0),V)==0).
					collect(Collectors.toList());
			if(L2.isEmpty()) throw new IllegalStateException("???");
			if(L2.size()==1) return L2.get(0);
			final VariantContext ctx = L2.stream().sorted((V1,V2)->{
				if(V1.isVariant() && !V2.isVariant()) return -1;
				if(!V1.isVariant() && V2.isVariant()) return 1;
				if(!V1.isFiltered() && V2.isFiltered()) return -1;
				if(V1.isFiltered() && !V2.isFiltered()) return 1;
				if(V1.isSNP() && !V2.isSNP()) return -1;
				if(!V1.isSNP() && V2.isSNP()) return 1;
				return 0;
				}).findFirst().get();
			if(debug_multiple<10)
				{
				LOG.warn("Multiple CHROM/POS/REF in "+this.uri+" at "+ctx.getContig()+":"+ctx.getStart()+":"+ctx.getReference().getDisplayString());
				debug_multiple++;
				}
			return ctx;
			}
		
		void reset(final VariantContext ctx0)
			{
			this.buffer.removeIf(V->
					 V.getContig().equals(ctx0.getContig()) &&
					 V.getStart() == ctx0.getStart()  &&
				     V.getReference().equals(ctx0.getReference()));
			}
		
		@Override
		public void close()
			{
			LOG.info("Closing "+this.uri+" Variants:"+this.count);
			CloserUtil.close(this.iter);
			}
		@Override
		public String toString() {
			return this.uri;
			}
		}
	
	
	
	@SuppressWarnings("resource")
	@Override
	public int doWork(final List<String> args) {
		if(this.archiveFile==null)
			{
			LOG.error("undefined output");
			return -1;
			}
		if(!this.archivePrefix.isEmpty() && !archivePrefix.endsWith("."))
			{
			this.archivePrefix=this.archivePrefix+".";
			}
		if(!this.archiveFile.getName().endsWith(".zip"))
			{
			IOUtil.assertDirectoryIsWritable(this.archiveFile);
			}
		
		
		
		PeekVCF vcfIterator1=null;
		PeekVCF vcfIterator2=null;
		ArchiveFactory archiveFactory=null;
		IntervalTreeMap<Boolean> capture = null;
		PrintWriter makefileWriter=null;
		try {
			if(args.size()==1)
				{
				LOG.info("Reading from VCF1=stdin and VCF2="+ args.get(0));
				vcfIterator1 = new PeekVCF(VCFUtils.createVCFIteratorFromInputStream(stdin()),"<STDIN>");
				vcfIterator2 = new PeekVCF(VCFUtils.createVCFIterator( args.get(0)),args.get(0));
				}
			else if(args.size()==2)
				{
				LOG.info("Reading from VCF1="+ args.get(0)+" and VCF2="+ args.get(1));
				vcfIterator1 = new PeekVCF(VCFUtils.createVCFIterator( args.get(0)), args.get(0));
				vcfIterator2 = new PeekVCF(VCFUtils.createVCFIterator( args.get(1)), args.get(1));
				}
			else
				{
				LOG.error("illegal number of arguments");
				return -1;
				}
				
			if( this.captureFile !=null )
				{
				LOG.info("Reading "+this.captureFile);
				capture = super.readBedFileAsBooleanIntervalTreeMap(this.captureFile);
				}
			
			this.global_dictionary = vcfIterator1.dict;
			if( !SequenceUtil.areSequenceDictionariesEqual( vcfIterator1.dict,  vcfIterator2.dict))
				{
				throw new JvarkitException.DictionariesAreNotTheSame( vcfIterator1.dict,  vcfIterator2.dict);
				}
			
			for(int side=0;side<2;++side)
				{
				final List<String> jexlExprStrings = (side==0?this.jexlExprStrings1:this.jexlExprStrings2);
				final PeekVCF peek = (side==0?vcfIterator1:vcfIterator2);
				//initialize JEXL map
				if(!jexlExprStrings.isEmpty()) continue;
					
				final ArrayList<String> dummyNames = new ArrayList<String>(jexlExprStrings.size());
			        for (int expCount =1; expCount< jexlExprStrings.size();++expCount ) {
			            dummyNames.add(String.format("vce%d",expCount));
			        	}
			    peek.jexlVCMatchExps.addAll(VariantContextUtils.initializeMatchExps(dummyNames, jexlExprStrings));
				}
			
			/* samples */
			final Set<String> samples0=new HashSet<>(vcfIterator1.header.getSampleNamesInOrder());
			final Set<String> samples1=new HashSet<>(vcfIterator2.header.getSampleNamesInOrder());
			final Set<String> commonSamples= new TreeSet<>(samples0);
			commonSamples.retainAll(samples1);
			
			if(commonSamples.size()!=samples0.size() || commonSamples.size()!=samples1.size())
				{
				LOG.warn("Warning: Not the same samples set. Using intersection of both lists.");
				}
			if(commonSamples.isEmpty())
				{	
				LOG.error("No common samples");
				return -1;
				}
			
			final Map<String,SampleInfo> sample2info=new HashMap<>(commonSamples.size());
			for(final String sampleName:commonSamples)
				{
				sample2info.put(sampleName, new  SampleInfo(sampleName));
				}
			
			
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(this.global_dictionary);
			for(;;)
				{
				VariantContext ctx0 = vcfIterator1.peek();
				VariantContext ctx1 = vcfIterator2.peek();
				
				final VariantContext smallest;
				if(ctx0==null && ctx1==null)
					{
					smallest=null;
					break;
					}
				else if(ctx0==null && ctx1!=null)
					{
					smallest = ctx1;
					}
				else if(ctx0!=null && ctx1==null)
					{
					smallest = ctx0;
					}
				else 
					{
					final int diff = this.compareChromPosRef.compare(ctx0,ctx1);
					if(diff<0)
						{
						smallest=ctx0;
						ctx1=null;
						}
					else if(diff>0)
						{
						smallest=ctx1;
						ctx0=null;
						}
					else
						{
						smallest=ctx0;
						}
					}
				progress.watch(smallest);
				vcfIterator1.reset(smallest);
				vcfIterator2.reset(smallest);
				
				if(capture!=null)
					{
					final Interval interval=  new Interval(smallest.getContig(),smallest.getStart(),smallest.getEnd());
					if(! capture.containsOverlapping(interval)) continue;
					}
				for(final String sampleName: sample2info.keySet())
					{
					sample2info.get(sampleName).visit(ctx0,ctx1);
					}
				}
			progress.finish();
			vcfIterator1.close();
			vcfIterator2.close();

		
			
			
			archiveFactory = ArchiveFactory.open(this.archiveFile);
			makefileWriter = archiveFactory.openWriter(this.archivePrefix+"Makefile");
			makefileWriter.println(".PHONY:all all2");
			makefileWriter.println("ALL_TARGETS=");
			makefileWriter.println("all:all2");
			
			for(final SampleInfo sampleInfo:sample2info.values())
				{
				for(final SampleCategory sampleCat:sampleInfo.sampleCat.values())
					{
					final String basename = this.archivePrefix + 
							sampleInfo.sampleName + "."+
							sampleCat.variantCatName.replace(' ', '_')
							;
					final String tsv =basename+".tsv";
					
					PrintWriter dataW = archiveFactory.openWriter(tsv);
					for(final String sk: new TreeSet<String>(sampleCat.counter.keySet()))
						{
						dataW.print(escapeUnderscore(sk));
						dataW.print("\t");
						dataW.print(sampleCat.counter.count(sk));
						dataW.println();
						}
					dataW.flush();dataW.close();
					
					
					final String png ="$(addsuffix .png,"+basename+")";
					makefileWriter.println("ALL_TARGETS+=" + png);
					
					makefileWriter.println(png+":"+tsv);
					makefileWriter.println("\techo '"
							+ "set ylabel \"Number of Genotypes  " + escapeUnderscore(sampleInfo.sampleName) +"\";"
							+ "set yrange [0:];"
							+ "set xlabel \"Category "+escapeUnderscore(vcf1Name)+": "+vcfIterator1.count+
									", "+escapeUnderscore(vcf2Name)+": "+vcfIterator2.count+" variants  \";"
							+ "set xtic rotate by 90 right;"
							+ "set size  ratio 0.618;"
							+ "set title \""+escapeUnderscore(vcf1Name)+" vs "+escapeUnderscore(vcf2Name)+" : Genotypes " + escapeUnderscore(sampleInfo.sampleName) +" / Variants: "+ escapeUnderscore(sampleCat.variantCatName) +" \";"
							+ "set key  off;"
							+ "set datafile separator \"\t\";"
							+ "set auto x;"
							+ "set style histogram;"
							+ "set style data histogram;"
							+ "set style fill solid border -1;"
							+ "set terminal png truecolor  size "+(500+sampleCat.counter.getCountCategories()*50)+", 1500;"
							+ "set output \"$@\";"
							+ "plot \"$<\" using 2:xtic(1) ti \"\";' | "
							+ "gnuplot");				
					}
				}
		

			makefileWriter.println("all2:${ALL_TARGETS}");

			makefileWriter.close();makefileWriter=null;
			archiveFactory.close();archiveFactory=null;
			return RETURN_OK;
			} 
		catch (final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(makefileWriter);
			CloserUtil.close(archiveFactory);
			CloserUtil.close(vcfIterator1);
			CloserUtil.close(vcfIterator1);
			}
		}

	
	
	public static void main(final String[] args) {
		new VcfCompareCallers().instanceMainWithExit(args);
	}
}
