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
package com.github.lindenb.jvarkit.tools.groupbygene;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.regex.Pattern;
import java.util.stream.Collectors;


import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.LineIterator;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

/**
BEGIN_DOC

## Motivation

Group VCF data by gene/transcript. By default it tries to use data from VEP and SnpEff

## Example

### Delimited output

```
$ curl -s -k "https://raw.github.com/arq5x/gemini/master/test/test4.vep.snpeff.vcf" |\
java -jar dist/groupbygene.jar |\
head | column  -t

#chrom  min.POS    max.POS    gene.name  gene.type         samples.affected  count.variations  M10475  M10478  M10500  M128215
chr10   52004315   52004315   ASAH2      snpeff-gene-name  2                 1                 0       0       1       1
chr10   52004315   52004315   ASAH2      vep-gene-name     2                 1                 0       0       1       1
chr10   52497529   52497529   ASAH2B     snpeff-gene-name  2                 1                 0       1       1       0
chr10   52497529   52497529   ASAH2B     vep-gene-name     2                 1                 0       1       1       0
chr10   48003992   48003992   ASAH2C     snpeff-gene-name  3                 1                 1       1       1       0
chr10   48003992   48003992   ASAH2C     vep-gene-name     3                 1                 1       1       1       0
chr10   126678092  126678092  CTBP2      snpeff-gene-name  1                 1                 0       0       0       1
chr10   126678092  126678092  CTBP2      vep-gene-name     1                 1                 0       0       0       1
chr10   135336656  135369532  CYP2E1     snpeff-gene-name  3                 2                 0       2       1       1
```

```

## History

* 201707: added pedigree, removed XML output


END_DOC
 */
@Program(
		name="groupbygene",
		keywords={"vcf","gene"},
		biostars={342790},
		description="Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff"
		)
public class GroupByGene
	extends Launcher
	{
	private static final Logger LOG = Logger.build(GroupByGene.class).make();

	@Parameter(names={"--filtered"},description="ignore FILTERED variants")
	private boolean ignore_filtered = false;
	@Parameter(names={"--gtFiltered"},description="[20170725] ignore FILTERED genotypes")
	private boolean ignore_filtered_genotype = false;
	@Parameter(names={"-T","--tag"},description="search Tag in INFO field containing the name of the genes (optional)")
	private Set<String> user_gene_tags = new HashSet<>();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outFile=null;
	@Parameter(names={"-p","--ped","--pedigree"},description="[20170725] "+ Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile=null;
	@Parameter(names={"--annIntergenic"},description="[20170726] Accept snpEff 'ANN' intergenic regions.")
	private boolean accept_snpeff_intergenic=false;
	@Parameter(names={"--typeRegexExclude"},description="[20170726] ignore prediction type matching this regex. e.g: '(ann_gene_id|ann_feature_id)' ")
	private String typeRegexExclude=null;
	@Parameter(names={"--slidingWindowSize"},description="[20170726] if greater than 0, add a sliding window of this size as the name of a virtual region ")
	private int slidingWindowSize=0;
	@Parameter(names={"--slidingWindowShift"},description="[20170726] if greater than 0, add shift the sliding window by this distance ")
	private int slidingWindowShift=0;
	@Parameter(names={"--fisher"},description="[20170726] Print fisher for case/control (experimental, need to work on this)")
	private boolean print_fisher=false;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
	
	/** the SAMSequenceDictionary used to sort reference */
	private SAMSequenceDictionary the_dictionary = null;
	/** the VCFCodec used to serialize variant */
	private AbstractVCFCodec the_codec = null;


	private final Function<String, Integer> contig2tid = (S)->{
		final int tid = the_dictionary.getSequenceIndex(S);
		if(tid<0) throw new JvarkitException.ContigNotFoundInDictionary(S, the_dictionary);
		return tid;
		};
	
	private final Comparator<String> contigComparator = (S1,S2) -> {
		if(S1.equals(S2)) return 0;
		if(the_dictionary==null || the_dictionary.isEmpty()) {
			return S1.compareTo(S2);
			} else {
				return contig2tid.apply(S1) - contig2tid.apply(S2);
			}
		};

	
	
	
	private static class GeneName
		{
		final String name;
		final String type;
		GeneName(final String name,final String type)
			{
			this.name=name;
			this.type=type;
			}
		@Override
		public int hashCode()
			{
			final int prime = 31;
			int result = 1;
			result = prime * result +  name.hashCode();
			result = prime * result +  type.hashCode();
			return result;
			}
		@Override
		public boolean equals(final Object o)
			{
			if (this == o) return true;
			if (o == null) return false;
			if (getClass() != o.getClass()) return false;
			final GeneName g=(GeneName)o;
			return name.equals(g.name) && type.equals(g.type);
			}
		@Override
		public String toString() {
			return  name+"("+type+")";
			}
		
		}
	
	private class Call implements Comparable<Call>
		{
		GeneName gene;
		String line;
		
		
		String getContig()
			{
			final int tab=this.line.indexOf(VCFConstants.FIELD_SEPARATOR);
			if(tab<1) throw new IllegalStateException("Cannot find tab in "+this.line);
			return this.line.substring(0, tab);
			}
		
		@Override
		public int compareTo(final Call o) {
			int i= contigComparator.compare(this.getContig(),o.getContig());
			if(i!=0) return i;
			i= this.gene.name.compareTo(o.gene.name);
			if(i!=0) return i;
			i= this.gene.type.compareTo(o.gene.type);
			return i;
			}
		}
	
	private class CallCodec
		extends AbstractDataCodec<Call>
		{
		@Override
		public void encode(final DataOutputStream dos,final Call c)
				throws IOException
			{
			dos.writeUTF(c.gene.name);
			dos.writeUTF(c.gene.type);
			writeString(dos, c.line);
			}
		
		@Override
		public Call decode(final DataInputStream dis) throws IOException
			{
			final String gName;
			try {
				gName=dis.readUTF();
			} catch (final Exception e) {
				return null;
				}
			final String gType=dis.readUTF();
			final Call c= new Call();
			c.gene=new GeneName(gName, gType);
			c.line = readString(dis);
			return c;
			}
		@Override
		public CallCodec clone() {
			return new CallCodec();
			}
		}
	
	
		
	
	
	public GroupByGene()
		{
		}
	
	private Set<GeneName> getGenes(final VcfTools vcfTools,final VariantContext ctx)
		{
		final Set<GeneName> set=new HashSet<GeneName>();
		for(final VepPredictionParser.VepPrediction pred: vcfTools.getVepPredictions(ctx))
			{
			String s=pred.getGeneName();
			if(s!=null)  set.add(new GeneName(s,"vep_gene_name"));
			s=pred.getEnsemblGene();
			if(s!=null)  set.add(new GeneName(s,"vep_ensembl_gene_name"));
			s=pred.getSymbol();
			if(s!=null)  set.add(new GeneName(s,"vep_symbol"));
			s=pred.getHGNC();
			if(s!=null)  set.add(new GeneName(s,"vep_hgnc"));
			s=pred.getHgncId();
			if(s!=null)  set.add(new GeneName(s,"vep_hgnc_id"));
			s=pred.getGene();
			if(s!=null)  set.add(new GeneName(s,"vep_gene"));
			s=pred.getRefSeq();
			if(s!=null)  set.add(new GeneName(s,"vep_refseq"));
			}
		for(final SnpEffPredictionParser.SnpEffPrediction pred: vcfTools.getSnpEffPredictions(ctx))
			{
			String s=pred.getGeneName();
			if(s!=null)  set.add(new GeneName(s,"snpeff_gene_name"));
			s=pred.getEnsemblTranscript();
			if(s!=null)  set.add(new GeneName(s,"snpeff_ensembl_transcript_name"));
			}
		for(final AnnPredictionParser.AnnPrediction pred: vcfTools.getAnnPredictions(ctx)) {
			if(pred.isIntergenicRegion() && !this.accept_snpeff_intergenic) continue;
			String s=pred.getGeneId();
			if(s!=null)  set.add(new GeneName(s,"ann_gene_id"));
			s=pred.getGeneName();
			if(s!=null)  set.add(new GeneName(s,"ann_gene_name"));
			s=pred.getFeatureId();
			if(s!=null)  set.add(new GeneName(s,"ann_feature_id"));
			}
		
		for(final String user_gene_tag:user_gene_tags)
			{
			if(user_gene_tag.isEmpty()) continue;
			if(user_gene_tag.equals(".")) continue;
			for(final String t:ctx.getAttributeAsList(user_gene_tag).stream().
					filter(O->O!=null).
					map(O->O.toString()).
					filter(S->!S.isEmpty()).
					collect(Collectors.toSet())
					)
				{
				set.add(new GeneName(t,"user:"+user_gene_tag));
				}
			}
		
		if(this.slidingWindowSize>0 && this.slidingWindowShift>0)
			{
			int p = ctx.getStart() - ctx.getStart()%this.slidingWindowSize;
			while(p-this.slidingWindowShift+this.slidingWindowSize >= ctx.getStart()  )
				{
				p-=this.slidingWindowShift;
				}
			p = Math.max(0 /* yes 0 */, p);
			while(p<=ctx.getStart() && p+this.slidingWindowSize >=ctx.getStart())
				{
				final String winName= ctx.getContig()+"_"+
						String.format("%06d",(p+1))+"_"+
						String.format("%06d",(p+this.slidingWindowSize));
				set.add(new GeneName(winName,"sliding_"+this.slidingWindowSize+"_"+this.slidingWindowShift));
				p+=this.slidingWindowShift;
				}
			
			}
		
		set.removeIf(G->G.name.isEmpty());
		return set;
		}
	

	private void read(final String input) throws IOException
		{
		LineIterator lineiter=null;
		SortingCollection<Call> sortingCollection=null;
		
		
		
		try {
			final Pattern regexType=(StringUtil.isBlank(this.typeRegexExclude)?null:Pattern.compile(this.typeRegexExclude));
			
			lineiter = (input==null?
						IOUtils.openStreamForLineIterator(stdin()):
						IOUtils.openURIForLineIterator(input)
						);

			sortingCollection =SortingCollection.newInstance(
					Call.class,
					new CallCodec(),
					(C1,C2)->{
						int i= C1.compareTo(C2);
						if(i!=0) return i;
						return C1.line.compareTo(C2.line);
					},
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sortingCollection.setDestructiveIteration(true);
	
			
			final VCFUtils.CodecAndHeader cah =VCFUtils.parseHeader(lineiter);
			final VCFHeader header = cah.header;
			this.the_dictionary = header.getSequenceDictionary();
			if(this.the_dictionary==null || this.the_dictionary.isEmpty())
				{
				throw new JvarkitException.DictionaryMissing(input);
				}
			this.the_codec = cah.codec;
			
			final List<String> sampleNames;
			if(header.getSampleNamesInOrder()!=null)
				{
				sampleNames = header.getSampleNamesInOrder();
				}
			else
				{
				sampleNames = Collections.emptyList();
				}
			
			final VcfTools vcfTools = new VcfTools(header);
			final Pedigree pedigree;
			if(this.pedigreeFile!=null) {
				pedigree = Pedigree.newParser().parse(this.pedigreeFile);
				}
			else
				{
				pedigree = Pedigree.newParser().parse(header);
				}
			
			final CharSplitter tab = CharSplitter.TAB;
			final ProgressFactory.Watcher<VariantContext> progress= ProgressFactory.newInstance().dictionary(the_dictionary).logger(LOG).build();
			while(lineiter.hasNext())
				{
				String line = lineiter.next();
				final VariantContext ctx = progress.apply(this.the_codec.decode(line));
				if(!ctx.isVariant()) continue;
				if(ignore_filtered && ctx.isFiltered()) continue;
				
				//simplify line
				final String tokens[]=tab.split(line);
				tokens[2]=VCFConstants.EMPTY_ID_FIELD;//ID
				tokens[5]=VCFConstants.MISSING_VALUE_v4;//QUAL
				tokens[6]=VCFConstants.UNFILTERED;//FILTER
				tokens[7]=VCFConstants.EMPTY_INFO_FIELD;//INFO
				line = String.join(VCFConstants.FIELD_SEPARATOR, Arrays.asList(tokens));
				
				for(final GeneName g:getGenes(vcfTools,ctx))
					{
					if(regexType!=null && regexType.matcher(g.type).matches()) continue;
					final Call c=new Call();
					c.line=line;
					c.gene=g;
					sortingCollection.add(c);
					}
				}
			CloserUtil.close(lineiter);lineiter=null;
			sortingCollection.doneAdding();
			progress.close();
			
			
			/** dump */			
	
			final Set<String> casesSamples = pedigree.getPersons().stream().
						filter(P->P.isAffected()).
						map(P->P.getId()).
						filter(ID->sampleNames.contains(ID)).
						collect(Collectors.toSet())
						;
			final Set<String> controlsSamples = pedigree.getPersons().stream().
					filter(P->P.isUnaffected()).
					map(P->P.getId()).
					filter(ID->sampleNames.contains(ID)).
					collect(Collectors.toSet())
					;
			
			final Set<String> maleSamples = pedigree.getPersons().stream().
					filter(P->P.isMale()).
					map(P->P.getId()).
					filter(ID->sampleNames.contains(ID)).
					collect(Collectors.toSet())
					;
			final Set<String> femaleSamples = pedigree.getPersons().stream().
					filter(P->P.isFemale()).
					map(P->P.getId()).
					filter(ID->sampleNames.contains(ID)).
					collect(Collectors.toSet())
					;
			final Predicate<Genotype> genotypeFilter = genotype -> {
				if(!genotype.isAvailable()) return false;
				if(!genotype.isCalled()) return false;
				if(genotype.isNoCall()) return false;
				if(genotype.isHomRef()) return false;
				if(this.ignore_filtered_genotype && genotype.isFiltered()) return false;
				return true;
			};
			
			
			
			PrintStream pw = openFileOrStdoutAsPrintStream(this.outFile);
			
			
			pw.print("#chrom");
			pw.print('\t');
			pw.print("min.POS");
			pw.print('\t');
			pw.print("max.POS");
			pw.print('\t');
			pw.print("gene.name");
			pw.print('\t');
			pw.print("gene.type");
			pw.print('\t');			
			pw.print("samples.affected");
			pw.print('\t');
			pw.print("count.variations");
			
			if(!casesSamples.isEmpty())
				{
				pw.print('\t');
				pw.print("pedigree.cases");
				}
			if(!controlsSamples.isEmpty())
				{
				pw.print('\t');
				pw.print("pedigree.controls");
				}
			if(!maleSamples.isEmpty())
				{
				pw.print('\t');
				pw.print("pedigree.males");
				}
			if(!femaleSamples.isEmpty())
				{
				pw.print('\t');
				pw.print("pedigree.females");
				}
			
			if(this.print_fisher && !controlsSamples.isEmpty() && !casesSamples.isEmpty())
				{
				pw.print('\t');
				pw.print("fisher");
				}

			
			for(final String sample:sampleNames)
				{
				pw.print('\t');
				pw.print(sample);
				}
			
			pw.println();
				
			
			final CloseableIterator<Call> iter=sortingCollection.iterator();
			final EqualRangeIterator<Call> eqiter = new EqualRangeIterator<>(iter, (C1,C2)->C1.compareTo(C2));
			while(eqiter.hasNext())
				{
				final List<Call> row = eqiter.next();
				final Call first= row.get(0);
	
				final List<VariantContext> variantList =row.stream().map(R-> GroupByGene.this.the_codec.decode(R.line)).collect(Collectors.toList());
				final int minPos = variantList.stream().mapToInt(R->R.getStart()).min().getAsInt();
				final int maxPos = variantList.stream().mapToInt(R->R.getEnd()).max().getAsInt();
				final Set<String> sampleCarryingMut = new HashSet<String>();
				final Counter<String> pedCasesCarryingMut = new Counter<String>();
				final Counter<String> pedCtrlsCarryingMut = new Counter<String>();
				final Counter<String> malesCarryingMut = new Counter<String>();
				final Counter<String> femalesCarryingMut = new Counter<String>();
				final Counter<String> sample2count = new Counter<String>();
				for(final VariantContext ctx: variantList)
					{
					for(final Genotype genotype:ctx.getGenotypes())
						{
						if(!genotypeFilter.test(genotype)) continue;
	 					final String sampleName= genotype.getSampleName();
	 					
						sample2count.incr(sampleName);
						sampleCarryingMut.add(sampleName);
						
						if(casesSamples.contains(sampleName))
							{
							pedCasesCarryingMut.incr(sampleName);
							}
						if(controlsSamples.contains(sampleName))
							{
							pedCtrlsCarryingMut.incr(sampleName);
							}
						if(maleSamples.contains(sampleName))
							{
							malesCarryingMut.incr(sampleName);
							}
						if(femaleSamples.contains(sampleName))
							{
							femalesCarryingMut.incr(sampleName);
							}
						}
					}
				
				
					pw.print(first.getContig());
					pw.print('\t');
					pw.print(minPos-1);//convert to bed
					pw.print('\t');
					pw.print(maxPos);
					pw.print('\t');
					pw.print(first.gene.name);
					pw.print('\t');
					pw.print(first.gene.type);
					pw.print('\t');
					pw.print(sampleCarryingMut.size());
					pw.print('\t');
					pw.print(variantList.size());
					
					if(!casesSamples.isEmpty())
						{
						pw.print('\t');
						pw.print(pedCasesCarryingMut.getCountCategories());
						}
					if(!controlsSamples.isEmpty())
						{
						pw.print('\t');
						pw.print(pedCtrlsCarryingMut.getCountCategories());
						}
					if(!maleSamples.isEmpty())
						{
						pw.print('\t');
						pw.print(malesCarryingMut.getCountCategories());
						}
					if(!femaleSamples.isEmpty())
						{
						pw.print('\t');
						pw.print(femalesCarryingMut.getCountCategories());
						}
					
					if(this.print_fisher && !controlsSamples.isEmpty() && !casesSamples.isEmpty())
						{
						int count_case_mut =0;
						int count_ctrl_mut = 0;
						int count_case_wild = 0;
						int count_ctrl_wild = 0;
						
						for(final VariantContext ctx: variantList) {
							for(final Genotype genotype:ctx.getGenotypes())
								{
								final String sampleName =  genotype.getSampleName();
								final boolean has_mutation = genotypeFilter.test(genotype);
								if( controlsSamples.contains(sampleName)) {
									if(has_mutation)
										{
										count_ctrl_mut++;
										}
									else
										{
										count_ctrl_wild++;
										}
									}
								else if( casesSamples.contains(sampleName)) {
									if(has_mutation)
										{
										count_case_mut++;
										}
									else
										{
										count_case_wild++;
										}
									}
								}
							}
						
						final FisherExactTest fisher = FisherExactTest.compute(
								count_case_mut,count_case_wild,
								count_ctrl_mut,count_ctrl_wild
								);
						pw.print('\t');
						pw.print(fisher.getAsDouble());
						}
	
					
					for(final String sample: sampleNames)
						{
						pw.print('\t');
						pw.print(sample2count.count(sample));
						}
					pw.println();
					if(pw.checkError()) break;
					
				}
			eqiter.close();
			iter.close();
			pw.flush();
			if(this.outFile!=null) pw.close();
			
			
			}
		finally
			{
			CloserUtil.close(lineiter);
			if(sortingCollection!=null) sortingCollection.cleanup();
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.slidingWindowSize>0)
			{
			if(this.slidingWindowShift>this.slidingWindowSize)
				{
				LOG.error("bad slidingw window parameters: shift >size ");
				return -1;
				}
			}
		
		try
			{
			read(oneFileOrNull(args));
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		}

	public static void main(final String[] args)
		{
		new GroupByGene().instanceMainWithExit(args);
		}
	}
