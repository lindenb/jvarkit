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
package com.github.lindenb.jvarkit.tools.groupbygene;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;


import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory;
import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;

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
		description="Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff",
		modificationDate="20210520",
		creationDate="20131209"
		)
public class GroupByGene
	extends Launcher
	{
	private static final Logger LOG = Logger.build(GroupByGene.class).make();

	@Parameter(names={"--filtered"},description="ignore FILTERED variants")
	private boolean ignore_filtered = false;
	@Parameter(names={"--gtFiltered"},description="[20170725] ignore FILTERED genotypes")
	private boolean ignore_filtered_genotype = false;
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile=null;
	@Parameter(names={"-p","--ped","--pedigree"},description="[20170725] "+ PedigreeParser.OPT_DESC)
	private Path pedigreePath=null;
	@Parameter(names={"-positions"},description="include variants positions in the output table.")
	private boolean print_positions=false;

	@Parameter(names={"--fisher"},description="[20170726] Print fisher for case/control (experimental, need to work on this)")
	private boolean print_fisher=false;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	@Parameter(names={"-l","--list"},description= "[20190626]list all available gene extractors", help=true)
	private boolean list_extractors = false;
	@Parameter(names={"-e","-E","--extractors"},description="[20190626]"+GeneExtractorFactory.OPT_DESC)
	private String extractorsNames="ANN/GeneId VEP/GeneId BCSQ/gene SMOOVE";

	
	/** the SAMSequenceDictionary used to sort reference */
	private ContigDictComparator contigDictComparator = null;



	
	
	
	private static class GeneName
		{
		final String name;
		final String label;
		final String type;
		GeneName(final String name,final String label,final String type)
			{
			this.name=name;
			this.label=StringUtils.isBlank(label)?".":label;
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
		VariantContext ctx;
		
		
		String getContig()
			{
			return ctx.getContig();
			}
		
		@Override
		public int compareTo(final Call o) {
			int i= contigDictComparator.compare(this.getContig(),o.getContig());
			if(i!=0) return i;
			i= this.gene.name.compareTo(o.gene.name);
			if(i!=0) return i;
			i= this.gene.type.compareTo(o.gene.type);
			return i;
			}
		
		public int compare2(final Call C2) {
			int i= this.compareTo(C2);
			if(i!=0) return i;
			i =  this.ctx.getContig().compareTo(C2.ctx.getContig());
			if(i!=0) return i;
			i =  Integer.compare(this.ctx.getStart(),C2.ctx.getStart());
			if(i!=0) return i;
			i =  this.ctx.getReference().compareTo(C2.ctx.getReference());
			return i;
			}
		
		
		}
	
	private class CallCodec 
		extends AbstractDataCodec<Call>
		{
		final VCFHeader header;
		private final VCFCodec vCodec;
		private final VCFEncoder vcfEncoder;

		CallCodec(final VCFHeader header) {
			this.header= header;
			this.vCodec = new VCFCodec();
			this.vcfEncoder = new VCFEncoder(header, false, false);
			this.vCodec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
			}
		@Override
		public void encode(final DataOutputStream dos,final Call c)
				throws IOException
			{
			dos.writeUTF(c.gene.name);
			dos.writeUTF(c.gene.label);
			dos.writeUTF(c.gene.type);
			writeString(dos, this.vcfEncoder.encode(c.ctx));
			}
		
		@Override
		public Call decode(final DataInputStream dis) throws IOException
			{
			final String gName;
			try {
				gName=dis.readUTF();
			} catch (final EOFException e) {
				return null;
				}
			final String gLbl=dis.readUTF();
			final String gType=dis.readUTF();
			final Call c= new Call();
			c.gene=new GeneName(gName,gLbl, gType);
			c.ctx = this.vCodec.decode(readString(dis));
			return c;
			}
		@Override
		public CallCodec clone() {
			return new CallCodec(this.header);
			}
		}
	
	
		
	
	
	public GroupByGene()
		{
		}
	
	private void read(final String input) throws IOException
		{
		SortingCollection<Call> sortingCollection=null;
		VCFIterator vcfIterator = null;
		
		
		try {
			final BcfIteratorBuilder iterbuilder = new BcfIteratorBuilder();
			vcfIterator = (input==null?
						iterbuilder.open(stdin()):
						iterbuilder.open(input)
						);
			final VCFHeader header = vcfIterator.getHeader();
			this.contigDictComparator = new ContigDictComparator(SequenceDictionaryUtils.extractRequired(header));


			sortingCollection =SortingCollection.newInstance(
					Call.class,
					new CallCodec(header),
					(C1,C2)->C1.compare2(C2),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sortingCollection.setDestructiveIteration(true);
	
			
			
			final GeneExtractorFactory geneExtractorFactory = new GeneExtractorFactory(header);
			final List<GeneExtractorFactory.GeneExtractor> geneExtractors = geneExtractorFactory.parse(this.extractorsNames);

			
			
			final List<String> sampleNames;
			if(header.getSampleNamesInOrder()!=null)
				{
				sampleNames = header.getSampleNamesInOrder();
				}
			else
				{
				sampleNames = Collections.emptyList();
				}
			
			final Pedigree pedigree;
			if(this.pedigreePath!=null) {
				final PedigreeParser pedParser = new PedigreeParser();
				pedigree = pedParser.parse(this.pedigreePath);
				}
			else
				{
				pedigree = PedigreeParser.empty();
				}
			
			final ProgressFactory.Watcher<VariantContext> progress= ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			while(vcfIterator.hasNext())
				{
				final VariantContext ctx = progress.apply(vcfIterator.next());
				if(!ctx.isVariant()) continue;
				if(ignore_filtered && ctx.isFiltered()) continue;
				
				//simplify line
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				vcb.noID();
				vcb.log10PError(VariantContext.NO_LOG10_PERROR);
				vcb.unfiltered();
				vcb.attributes(Collections.emptyMap());
				final VariantContext ctx2 = vcb.make();
				
				final SortingCollection<Call> finalSorter = sortingCollection; 
				geneExtractors.stream().
					flatMap(EX->EX.apply(ctx).keySet().stream()).
					forEach(KG->{
						final Call c=new Call();
						c.ctx=ctx2;
						c.gene=new GeneName(KG.getKey(),KG.getGene(),KG.getMethod());
						finalSorter.add(c);
					});
				
				}
			CloserUtil.close(vcfIterator);vcfIterator=null;
			sortingCollection.doneAdding();
			progress.close();
			
			
			/** dump */			
	
			final Set<String> casesSamples = pedigree.getAffectedSamples().stream().
						map(P->P.getId()).
						filter(ID->sampleNames.contains(ID)).
						collect(Collectors.toSet())
						;
			final Set<String> controlsSamples = pedigree.getUnaffectedSamples().stream().
					map(P->P.getId()).
					filter(ID->sampleNames.contains(ID)).
					collect(Collectors.toSet())
					;

			final Set<String> maleSamples = pedigree.getSamples().stream().
					filter(P->P.isMale()).
					map(P->P.getId()).
					filter(ID->sampleNames.contains(ID)).
					collect(Collectors.toSet())
					;
			final Set<String> femaleSamples = pedigree.getSamples().stream().
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
			
			
			
			PrintStream pw = openPathOrStdoutAsPrintStream(this.outFile);
			
			
			pw.print("#chrom");
			pw.print('\t');
			pw.print("min.POS");
			pw.print('\t');
			pw.print("max.POS");
			pw.print('\t');
			pw.print("gene.name");
			pw.print('\t');
			pw.print("gene.label");
			pw.print('\t');
			pw.print("gene.type");
			pw.print('\t');			
			pw.print("samples.affected");
			pw.print('\t');
			pw.print("count.variations");
			if(this.print_positions) {
				pw.print('\t');
				pw.print("positions");
				}
			
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
	
				final List<VariantContext> variantList =row.stream().map(R-> R.ctx).collect(Collectors.toList());
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
					pw.print(first.gene.label);
					pw.print('\t');
					pw.print(first.gene.type);
					pw.print('\t');
					pw.print(sampleCarryingMut.size());
					pw.print('\t');
					pw.print(variantList.size());
					if(this.print_positions) {
						pw.print('\t');
						pw.print(variantList.stream().
								map(CTX->String.valueOf(CTX.getStart())+":"+CTX.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/"))).
								collect(Collectors.joining(";"))
							);
						}
					
					
					
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
						
						for(final String sampleName : header.getSampleNamesInOrder())  {
							final boolean has_mutation  = variantList.stream().
									map(V->V.getGenotype(sampleName)).
									anyMatch(G->G!=null && genotypeFilter.test(G));
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
			CloserUtil.close(vcfIterator);
			if(sortingCollection!=null) sortingCollection.cleanup();
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
