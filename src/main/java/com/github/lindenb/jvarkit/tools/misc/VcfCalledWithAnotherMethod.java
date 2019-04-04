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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.PeekIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
/**
BEGIN_DOC

## Example


```make
SHELL=/bin/bash

define ff
dir1/dir2/sample_variations.$(1).annotations.vcf.gz
endef

all :
	java -jar  dist/vcfcalledwithanothermethod.jar \
		-f $(call ff,samtools) \
		-f $(call ff,varscan) \
		-f $(call ff,freebayes) \
			$(call ff,gatkHapCaller)
	
```

output:

```
(...)
1	12718	.	G	C	197.77	VariantNotFoundElseWhere	AC=1;AF=0.500;AN=2;BaseQRankSum=-1.418;ClippingRankSum=1.220;DP=22;QD=8.99;ReadPosRankSum=1.022;SEGDUP=1:10485-19844,1:10464-40733,1:10000-19844,1:10485-40733,1:10000-87112,1:10000-20818	GT:AD:COUNT_DISCORDANT:COUNT_SAME:DP:GQ:PL	0/1:12,10:0:0:22:99:226,0,286
1	23119	.	T	G	637.77	PASS	FOUND_COUNT=2;FOUND_KEY=sample_variations.varscan.annotations,sample_variations.samtools.annotations;FS=34.631;GERP_SCORE=-0.558;MLEAC=1;MLEAF=0.500;MQ=25.98;MQ0=0;MQRankSum=-2.888;POLYX=1;PRED=uc010nxq.1|||||intron_variant;QD=18.22;ReadPosRankSum=1.634;SEGDUP=1:10485-19844,1:10464-40733,1:10000-19844,1:10485-40733,1:10000-87112,1:10000-20818	GT:AD:COUNT_DISCORDANT:COUNT_SAME:DP:GQ:PL	0/1:17,18:0:2:35:99:666,0,727

```


END_DOC
 */
@Program(name="vcfcalledwithanothermethod",
		description="Compare one vcf with other , add a flag to tell if a variant was called with another method. Vcf must be sorted on the same Dict.",
		keywords={"vcf","compare","concordance"}
		)
public class VcfCalledWithAnotherMethod extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfCalledWithAnotherMethod.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-f","--vcfs"},description="Add alternate VCF files. File ending with '.list' will be interpreted as a list of path of vcf.")
	private Set<String> externalVcfStrs = new HashSet<>();
	
	@Parameter(names={"--filter"},description="FILTER name: the variant was NOT found in another VCF (CONTIG/POS/REF/at-least-one-ALT). Empty: no filter")
	private String filterNameVariantNotFoundElseWhere="VariantNotFoundElseWhere";
	
	@Parameter(names={"--foundKey"},description="INFO name for the file identifiers where a variant was found")
	private String infoFoundKey="FOUND_KEY";
	@Parameter(names={"--foundCount"},description="INFO name for the file identifiers where a variant was found")
	private String infoFoundCount="FOUND_COUNT";
	@Parameter(names={"--gtDiscordant"},description="FORMAT name for the number of time we didn't find the same genotype")
	private String formatCountDiscordant="COUNT_DISCORDANT";
	@Parameter(names={"--gtSame"},description="FORMAT name for the number of time we found the same genotype")
	private String formatCountSame="COUNT_SAME";
	@Parameter(names={"--nocallhomref"},description="NO_CALL is same as HOM_REF")
	private boolean noCallSameAsHomRef=false;

	
	
	private SAMSequenceDictionary dictionary=null;
	
	
	private final Function<String,Integer> contig2tid = C->{
		int i = dictionary.getSequenceIndex(C);
		if(i==-1) throw new JvarkitException.ContigNotFoundInDictionary(C,dictionary);
		return i;
		};
	private final Comparator<String> compareContig = (C1,C2)-> contig2tid.apply(C1)-contig2tid.apply(C2);
	
	
	private final Comparator<VariantContext> compareContextChromPos = (VC1,VC2)-> {
		int i = compareContig.compare(VC1.getContig(), VC2.getContig());
		if(i!=0) return i;
		return VC1.getStart() - VC2.getStart();
		};
	
	private final Comparator<VariantContext> compareContextChromPosRef = (VC1,VC2)-> {
		int i = compareContextChromPos.compare(VC1,VC2);
		if(i!=0) return i;
		return VC1.getReference().compareTo(VC2.getReference());
		};
	
	
	private class ExternalVcf
		{
		private final File vcfFile;
		private VCFFileReader reader;
		private VariantContext prev_ctx=null;
		private final List<VariantContext> buffer=new ArrayList<>();
		final CloseableIterator<VariantContext> citer;
		final PeekIterator<VariantContext> iter;
		String key;//unique key to identify this file in the VCF
		ExternalVcf(final File vcfFile)
			{
			this.vcfFile = vcfFile;
			this.reader = new VCFFileReader(this.vcfFile, false);
			this.citer = this.reader.iterator();
			this.iter = new PeekIterator<>(this.citer);
			this.key = VCFUtils.escapeInfoField(this.vcfFile.getName());
			if(key.endsWith(".gz")) key = key.substring(0, key.length()-3);
			if(key.endsWith(".vcf")) key = key.substring(0, key.length()-4);
			if(key.endsWith(".bcf")) key = key.substring(0, key.length()-4);
			
			if(this.reader.getFileHeader().getSequenceDictionary()==null) {
				throw new JvarkitException.VcfDictionaryMissing(this.vcfFile);
				}
			if(!SequenceUtil.areSequenceDictionariesEqual(this.reader.getFileHeader().getSequenceDictionary(), dictionary))
				{
				throw new JvarkitException.UserError("not same sequence dictionary between input and "+this.vcfFile);
				}
			}
		
		public List<VariantContext> _filter(List<VariantContext> L,final VariantContext ctx)
			{
			return L.isEmpty()?
				L:
				L.stream().
					filter(V->compareContextChromPosRef.compare(V, ctx)==0).
					collect(Collectors.toList());
			}			
		
		public List<VariantContext> get(final VariantContext ctx)
			{
			if(!this.buffer.isEmpty() &&  compareContextChromPos.compare(this.buffer.get(0), ctx)==0)
				{
				return _filter(this.buffer,ctx);
				}
			this.buffer.clear();
			for(;;)
				{
				if(!this.iter.hasNext())
					{
					break;
					}
				final VariantContext ctx2= this.iter.peek();
				if(this.prev_ctx!=null && compareContextChromPos.compare(this.prev_ctx, ctx2)>0)
					{
					throw new RuntimeIOException(
							"Bad order in "+ this.vcfFile+": got\n "+this.prev_ctx+"\nbefore\n "+ctx2+"\n"
							);
					}
				this.prev_ctx =  ctx2;
				final int d = compareContextChromPos.compare(ctx2, ctx);
				if( d < 0)
					{
					this.iter.next();//consumme
					continue;
					}
				else if(d>0)
					{
					break;
					}
				else //equal
					{
					this.buffer.add(this.iter.next());
					}
				}
			return _filter(this.buffer,ctx);
			}
		
		void close() {
			this.reader.close();
			this.citer.close();
			this.buffer.clear();
			}
		@Override
		public String toString() {
			return vcfFile.toString();
			}
		}
	
	private class GenotypeCount
		{	
		final Genotype gt;
		int countSame=0;
		int countDiscordant=0;
		
		GenotypeCount(final Genotype gt) { this.gt=gt;}
		
		}



	public VcfCalledWithAnotherMethod()
		{
		
		}
	@Override
	public int doWork(final List<String> args) {
		if(this.externalVcfStrs==null) {
			LOG.warn("Missing external VCFs");
			}
		return doVcfToVcf(args, outputFile);
		}
	
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator in,
			final VariantContextWriter out
			) {
		final List<ExternalVcf> externalVcfs = new ArrayList<>();
		try {
			final VCFHeader header= in.getHeader();
			this.dictionary = header.getSequenceDictionary();
			/** open primitive input */
			if(this.dictionary==null) {
				LOG.error("no dictionary in input");
				return -1;
				}
			
			final Set<File> samtoolsFiles = new HashSet<>();
			this.externalVcfStrs.stream().
					filter(S->S.endsWith(".list")).
					map(S->Paths.get(S)).
					forEach(P->{
						try 
						{
							samtoolsFiles.addAll(Files.readAllLines(P).stream().
								filter(L->!(L.startsWith("#") || L.trim().isEmpty())).
								map(S->new File(S)).
								collect(Collectors.toSet()));
						}
						catch(final Exception err)
						{
							throw new RuntimeIOException(err);
						}
					});
			
			samtoolsFiles.addAll(this.externalVcfStrs.stream().
					filter(S->!S.endsWith(".list")).
					map(S->new File(S)).
					collect(Collectors.toSet()));
			
			
			externalVcfs.addAll(samtoolsFiles.stream().
					map(F->new ExternalVcf(F)).
					collect(Collectors.toList())
					);
			/** check for uniq keys */
			final Set<String> uniqKeys=new HashSet<>();
			for(final ExternalVcf extvcf: externalVcfs)
				{
				int n=0;
				for(;;)
					{
					final String newkey=extvcf.key+(n==0?"":"_"+n);
					if(!uniqKeys.contains(newkey)) 
						{
						extvcf.key = newkey;
						uniqKeys.add(newkey);
						break;
						}
					++n;
					}
				}
			
			
			final VCFHeader h2=new VCFHeader(header);
			for(final ExternalVcf extvcf: externalVcfs)
				{
				h2.addMetaDataLine(new VCFHeaderLine("otherVcf:"+extvcf.key, extvcf.vcfFile.getPath()));
				}
			
			final VCFFilterHeaderLine variantNotFoundElsewhereFILTER = 
					(
					filterNameVariantNotFoundElseWhere.isEmpty() ?
					null:
					new VCFFilterHeaderLine(
						filterNameVariantNotFoundElseWhere,
						"Variant Was not found in other VCFs: "+externalVcfs.stream().map(S->S.vcfFile.getPath()).collect(Collectors.joining(", "))
							)
					);
			if(variantNotFoundElsewhereFILTER!=null)
				{
				h2.addMetaDataLine(variantNotFoundElsewhereFILTER);
				}
			
			final VCFInfoHeaderLine variantFoundElseWhereKeys = new VCFInfoHeaderLine(
					this.infoFoundKey,
					VCFHeaderLineCount.UNBOUNDED,
					VCFHeaderLineType.String,
					"Variant was found in the VCFs designed by those keys"
					);
			h2.addMetaDataLine(variantFoundElseWhereKeys);

			final VCFInfoHeaderLine variantFoundElseWhereCount = new VCFInfoHeaderLine(
					this.infoFoundKey,
					1,
					VCFHeaderLineType.Integer,
					"Number of times the Variant was found  in the VCFs"
					);
			h2.addMetaDataLine(variantFoundElseWhereCount);
			
			
			final VCFFormatHeaderLine genotypeCountSame = new VCFFormatHeaderLine(
					this.formatCountSame,
					1,
					VCFHeaderLineType.Integer,
					"Number of times the Genotype was found the same in other VCFs"
					);
			h2.addMetaDataLine(genotypeCountSame);
			final VCFFormatHeaderLine genotypeCountDiscordant = new VCFFormatHeaderLine(
					this.formatCountDiscordant,
					1,
					VCFHeaderLineType.Integer,
					"Number of times the Genotype was found dicordant in other VCFs"
					);
			h2.addMetaDataLine(genotypeCountDiscordant);
			
			super.addMetaData(h2);
			
			
			final VariantContextWriter w = super.openVariantContextWriter(outputFile);
			w.writeHeader(h2);
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(this.dictionary);
			while(in.hasNext())
				{
				final VariantContext ctx = progress.watch(in.next());
				final List<GenotypeCount> genotypeCounts = ctx.getGenotypes().
						stream().
						map(G->new GenotypeCount(G)).collect(Collectors.toList());
				
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				final Set<String> variantFoundInOtherVcfs = new HashSet<>();
				for(final ExternalVcf extvcf: externalVcfs)
					{
					final List<VariantContext> otherVariants = extvcf.get(ctx);
					
					if(otherVariants.stream().filter(CTX2->{
						if(ctx.getAlternateAlleles().isEmpty()) return true;
						for(final Allele a: ctx.getAlternateAlleles())
							{
							if(CTX2.hasAllele(a)) return true;
							}
						return false;
						}).findAny().isPresent())
						{
						variantFoundInOtherVcfs.add(extvcf.key);
						}
					for(final GenotypeCount gt: genotypeCounts)
						{
						for(final VariantContext  otherVariant : otherVariants)
							{
							final Genotype otherGt = otherVariant.getGenotype(gt.gt.getSampleName());
							if(otherGt==null) continue;
							if(gt.gt.sameGenotype(otherGt) || (this.noCallSameAsHomRef && ((gt.gt.isNoCall() && otherGt.isHomRef())|| (gt.gt.isHomRef() && otherGt.isNoCall())))) {
								gt.countSame++;
								}
							else
								{
								gt.countDiscordant++;
								}
							}
						}
					}
				
				
				vcb.genotypes(genotypeCounts.stream().map(G->{
					final GenotypeBuilder gb=new GenotypeBuilder(G.gt);
					gb.attribute(genotypeCountSame.getID(), G.countSame);
					gb.attribute(genotypeCountDiscordant.getID(), G.countDiscordant);
					return gb.make();
				}).collect(Collectors.toList()));
				
				
				vcb.attribute(
						variantFoundElseWhereCount.getID(),
						variantFoundInOtherVcfs.size()
						);
				if(variantFoundInOtherVcfs.isEmpty())
					{
					if(variantNotFoundElsewhereFILTER!=null)
						{
						vcb.filter(variantNotFoundElsewhereFILTER.getID());
						}
					}
				else
					{
					
					if(variantNotFoundElsewhereFILTER!=null && !ctx.isFiltered())
						{
						vcb.passFilters();
						}
					vcb.attribute(
							variantFoundElseWhereKeys.getID(),
							new ArrayList<>(variantFoundInOtherVcfs)
							);
					}
				
				w.add(vcb.make());
				}
			
			progress.finish();
			while(!externalVcfs.isEmpty()) externalVcfs.remove(0).close();

			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			while(!externalVcfs.isEmpty()) externalVcfs.remove(0).close();
			}
		}
	
	public static void main(String[] args) {
		new VcfCalledWithAnotherMethod().instanceMainWithExit(args);

	}

}
