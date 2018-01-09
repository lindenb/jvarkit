/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.lumpysv;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC

## Example

```
 java -jar dist/lumpyvcf2circos.jar  --minsu 50 -inv -bnb -dup  -o tmp  LumpyExpress.vcf.gz \
  && (cd tmp; /path/to/bin/circos  -outputdir . -conf lumpy.circos.conf  )
```

END_DOC
*/
@Program(name="lumpymoresamples",
		description="Lumpy to Circos",
		keywords={"lumpy","circos","sv","vcf"},
		generate_doc=false
		)
public class LumpyMoreSamples extends Launcher {
	
	private static final Logger LOG = Logger.build(LumpyMoreSamples.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-b","--bam"},description="File containing list of bam, one path per line.",required=true)
	private File bamFileList = null;
	
	private final Map<String,List<SamReader>> sample2samreaders = new HashMap<>();
	
	@Override
	public int doWork(final List<String> args) {
	VcfIterator r=null;
	VariantContextWriter vcw=null;
	try {
		r = super.openVcfIterator(oneFileOrNull(args));
		
		final VCFHeader headerIn =r.getHeader();
		final SAMSequenceDictionary dict = headerIn.getSequenceDictionary();
		if(dict==null)
			{
			LOG.error(JvarkitException.VcfDictionaryMissing.getMessage("input vcf"));
			return -1;
			}
		
		final SamReaderFactory samReaderFactory = SamReaderFactory.
				makeDefault().
				validationStringency(ValidationStringency.LENIENT);
		IOUtil.slurpLines(bamFileList).stream().forEach(F->{
			if(F.trim().isEmpty()) return;
			final SamReader sr = samReaderFactory.open(SamInputResource.of(F));
			final SAMFileHeader samHeader = sr.getFileHeader();
			final SAMSequenceDictionary dict2 = samHeader.getSequenceDictionary();
			if(dict2==null)
				{
				throw new JvarkitException.BamDictionaryMissing(F);
				}
			if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict2)) {
				throw new JvarkitException.DictionariesAreNotTheSame(dict,dict2);
				}
			for(final SAMReadGroupRecord rg: samHeader.getReadGroups())
				{
				final String sample = rg.getSample();
				List<SamReader> readers = this.sample2samreaders.get(sample);
				if(readers==null) {
					readers=new ArrayList<>();
					this.sample2samreaders.put(sample,readers);
					}
				readers.add(sr);
				}
			});
		
		final Set<String >genotypeSampleNames = new HashSet<>(headerIn.getSampleNamesInOrder());
		genotypeSampleNames.addAll(this.sample2samreaders.keySet());
		
		final VCFHeader headerOut =
				new VCFHeader(
						headerIn.getMetaDataInInputOrder(),
						genotypeSampleNames
						);

		vcw = super.openVariantContextWriter(this.outputFile);
		vcw.writeHeader(headerOut);
		while(r.hasNext())
			{
			final VariantContext ctx = r.next();
			final List<Genotype> genotypes = new ArrayList<>(ctx.getGenotypes());
			
			for(final String sample:this.sample2samreaders.keySet())
				{
				
				}
			
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			vcb.genotypes(genotypes);
			
			vcw.add(vcb.make());
			}
		r.close();r=null;
		this.sample2samreaders.values().stream().
			flatMap(L->L.stream()).
			forEach(R->CloserUtil.close(R));
		
		LOG.info("done");
		return 0;
		} 
	catch(final Exception err)
			{
				LOG.error(err);
				return -1;
			}
		finally {
			CloserUtil.close(r);		
		}
	}
	
	
public static void main(String[] args) {
	new LumpyMoreSamples().instanceMainWithExit(args);
}
}
