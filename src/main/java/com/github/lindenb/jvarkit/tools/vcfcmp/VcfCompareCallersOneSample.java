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
package com.github.lindenb.jvarkit.tools.vcfcmp;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC

## Misc.

I used 
```
java -jar dist/VcfCompareCallersOneSample.jar  -m1 -M1 -f samtools.vcf gatk.vcf
java -jar dist/VcfCompareCallersOneSample.jar  -m1 -M1 -f gatk.vcf samtools.vcf
```

shouldn't I get the same number of variants in both files ?
**Answer** is "not always" in the following case:

in gatk.vcf:

```
11	244197	rs1128322	T	C
````

in samtools.vcf:
```
11	244197	rs1128322	T	C,G
```

in
```
java -jar dist/VcfCompareCallersOneSample.jar  -m1 -M1 -f samtools.vcf gatk.vcf
```
we keep the variant because we found 'C' in samtools and 'gatk'

```
java -jar dist/VcfCompareCallersOneSample.jar  -m1 -M1 -f gatk.vcf samtools.vcf
```

the variant is discarded because 'G' is found in samtools but not in 'gatk.vcf'

END_DOC

 */
@Program(name="vcfcomparecallersonesample",
description="For my colleague Julien: VCF with one sample called using different callers. Only keep variant if it was found in min<x=other-files<=max",
keywords= {"vcf","compare"}
)
public class VcfCompareCallersOneSample
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfCompareCallersOneSample.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names="-f",description="VCF to be challenged.  Must be sorted on dict. Must contain a dict.")
	private Set<File> challengerVcf=new HashSet<File>();
	@Parameter(names="-m",description="min number of challengers found, inclusive.")
	private int minCountInclusive=0;
	@Parameter(names="-M",description=" max number of challengers found, inclusive.")
	private int maxCountInclusive=Integer.MAX_VALUE-1;
	@Parameter(names="-a",description="ignore ALT allele")
	private boolean ignoreAlternate=false;
	
	public VcfCompareCallersOneSample()
		{
		}
	
	public Set<File> getChallengerVcf() {
		return challengerVcf;
		}
	
	
	public void setMinCountInclusive(int minCountInclusive) {
		this.minCountInclusive = minCountInclusive;
		}
		
	public void setMaxCountInclusive(int maxCountInclusive) {
		this.maxCountInclusive = maxCountInclusive;
		}
	
	public void setIgnoreAlternate(boolean ignoreAlternate) {
		this.ignoreAlternate = ignoreAlternate;
		}
	
	
	
@Override
	public int doWork(final List<String> args) {
		File inputFile=null;
		final List<EqualRangeVcfIterator> listChallengers = new ArrayList<>();
		VariantContextWriter vcw=null;
		VCFIterator in=null;
		try {
			 in  = super.openVCFIterator(oneFileOrNull(args));
			
			
			VCFHeader header=in.getHeader();
			if(header.getNGenotypeSamples()!=1)
				{
				LOG.error("vcf must have only one sample");
				return -1;
				}
			
			final VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			h2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), h2);
			
			
			SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict==null)
				{
				LOG.error(JvarkitException.VcfDictionaryMissing.getMessage("input"));
				return -1;
				}
			final Comparator<VariantContext> ctxComparator = VCFUtils.createTidPosComparator(dict);

			/* load files to be challenged */
			for(final File cf :this.challengerVcf)
				{
				//do not challenge vs itself
				if(inputFile!=null && inputFile.equals(cf))
					{
					LOG.error("Ignoring challenger (self): "+cf);
					continue;
					}
				VCFIterator cin = VCFUtils.createVCFIteratorFromFile(cf);
				VCFHeader ch=cin.getHeader();
				if(ch.getNGenotypeSamples()!=1)
					{
					LOG.warning("vcf.must.have.only.one.sample");
					cin.close();
					continue;
					}
				if(!header.getSampleNamesInOrder().get(0).equals(
						ch.getSampleNamesInOrder().get(0)))
					{
					LOG.warning("Ignoring "+cf+" because not the same sample.");
					cin.close();
					continue;
					}
				SAMSequenceDictionary hdict = ch.getSequenceDictionary();
				if(hdict==null ||
					!SequenceUtil.areSequenceDictionariesEqual(dict, hdict))
					{
					LOG.error("not.the.same.sequence.dictionaries");
					return -1;
					}
				listChallengers.add(new EqualRangeVcfIterator(cin,ctxComparator));
				}
			
			vcw= super.openVariantContextWriter(outputFile);
			vcw.writeHeader(h2);
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
			VariantContext prev_ctx=null;
			while(in.hasNext() && !vcw.checkError())
				{
				VariantContext ctx = progress.watch(in.next());
				
				//check input order
				if(prev_ctx!=null && ctxComparator.compare(prev_ctx,ctx)>0)
					{
					LOG.error("bad sort order : got\n\t"+
							prev_ctx+"\nbefore\n\t"+
							ctx+"\n");
					return -1;
					}
				prev_ctx=ctx;
				
				
				int countInOtherFiles=0;
				for(final EqualRangeVcfIterator citer:listChallengers)
					{
					boolean foundInThatFile=false;
					final List<VariantContext> ctxChallenging = citer.next(ctx);
					for(final VariantContext ctx2:ctxChallenging)
						{
						if(!ctx2.getReference().equals(ctx.getReference())) continue;
						boolean ok=true;
						if(!this.ignoreAlternate)
							{
							Set<Allele> myAlt=new HashSet<Allele>(ctx.getAlternateAlleles());
							myAlt.removeAll(ctx2.getAlternateAlleles());
							if(!myAlt.isEmpty()) ok=false;
							}
						
						if(ok)
							{
							foundInThatFile=true;
							break;
							}
						}
					countInOtherFiles+=(foundInThatFile?1:0);
					}
				
				if(countInOtherFiles >= minCountInclusive &&
					countInOtherFiles <= maxCountInclusive)
					{
					vcw.add(ctx);
					}
				
				}
			progress.finish();
			return 0;
			} 
		catch (final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcw);
			CloserUtil.close(listChallengers);
			CloserUtil.close(in);
			}
		
		}
	
	
	public static void main(final String[] args) {
		new VcfCompareCallersOneSample().instanceMainWithExit(args);
	}
	}
