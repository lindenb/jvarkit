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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintStream;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example

```
bash
$   curl -sL "https://raw.githubusercontent.com/arq5x/bedtools2/bc2f97d565c36a82c1a0b12f570fed4398001e5f/test/map/test.vcf" |\
    java -jar dist/biostar130456.jar -x -z -p "sample.__SAMPLE__.vcf.gz" 
sample.NA00003.vcf.gz
sample.NA00001.vcf.gz
sample.NA00002.vcf.gz

$ gunzip -c sample.NA00003.vcf.gz
(...)
##source=myImputationProgramV3.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00003
chr1	10	rs6054257	G	A	29	PASS	AF=0.5;DB;DP=14;H2;NS=3	GT:DP:GQ:HQ	1/1:5:43
chr1	20	rs6040355	A	G,T	67	PASS	AA=T;AF=0.333,0.667;DB;DP=10;NS=2	GT:DP:GQ	2/2:4:35
chr1	130	microsat1	GTC	G,GTCT	50	PASS	AA=G;DP=9;NS=3	GT:DP:GQ	1/1:3:40
chr2	130	microsat1	GTC	G,GTCT	50	PASS	AA=G;DP=9;NS=3	GT:DP:GQ	1/1:3:40
```

## See also

 * GATK SelectVariants with option -sn , or bcftools view --samples-file

END_DOC
*/
@Program(
		name="biostar130456",
		description="Split individual VCF files from multisamples VCF file",
		keywords={"vcf","samples","sample"},
		biostars=130456
		)
public class Biostar130456 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar130456.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	private final static String SAMPLE_TAG="__SAMPLE__";
	@Parameter(names={"-p","--pattern"},description="output file pattern. Must contain the word "+SAMPLE_TAG,required=true)
	private String filepattern = null;

	@Parameter(names={"-x","--uncalled"},description="remove uncalled genotypes")
	private boolean remove_uncalled = false;
	@Parameter(names={"-z","--homref"},description="remove homzygote REF/REF")
	private boolean remove_homref = false;
	@Parameter(names={"-f","--filtered"},description="remove filtered Genotype")
	private boolean remove_filtered = false;
	
	@ParametersDelegate
	private VariantAttributesRecalculator recalculator = new VariantAttributesRecalculator();

	
	
	@Override
	public int doWork(final List<String> args) {
			if(StringUtils.isBlank(this.filepattern) || !filepattern.contains(SAMPLE_TAG))
				{
				LOG.error("File pattern is missing "+SAMPLE_TAG);
				return -1;
				}
			PrintStream out = null;
			VCFIterator in=null;
			final String inputName= oneFileOrNull(args);

			try
				{
				out = openFileOrStdoutAsPrintStream(outputFile);
				in = super.openVCFIterator(inputName);
				final VCFHeader header=in.getHeader();
				this.recalculator.setHeader(header);
				
				final Set<String> samples = new HashSet<String>(header.getSampleNamesInOrder());
				final Map<String,VariantContextWriter> sample2writer=new HashMap<>(samples.size());
	
				if(samples.isEmpty())
					{
					LOG.error("VCF doesn't contain any sample");
					return -1;
					}
				LOG.info("N sample:"+samples.size());
				for(final String sample:samples)
					{
					final VCFHeader h2=new VCFHeader(
							header.getMetaDataInInputOrder(),
							Collections.singleton(sample)
							);
					super.addMetaData(h2);
					final String sampleFile= filepattern.replaceAll(SAMPLE_TAG,sample);
					out.println(sampleFile);
					final File fout = new File(sampleFile);
					if(fout.getParentFile()!=null) fout.getParentFile().mkdirs();
					final VariantContextWriter w= VCFUtils.createVariantContextWriter(fout);
					w.writeHeader(h2);
					
					sample2writer.put(sample, w);
					}
				final ProgressFactory.Watcher<VariantContext> progress=ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
				while(in.hasNext())
					{
					final VariantContext ctx= progress.apply(in.next());
					for(final String sample: samples)
						{
						final Genotype g= ctx.getGenotype(sample);
						if(g==null) continue;
						if(remove_uncalled && !g.isCalled() ) continue;
						if(remove_homref && g.isHomRef()) continue;
						if(remove_filtered && g.isFiltered()) continue;
						final VariantContextWriter w= sample2writer.get(sample);
						final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
						final GenotypeBuilder gb=new GenotypeBuilder(g);
						final Set<Allele> alleles=new HashSet<>(3);
						alleles.add(ctx.getReference());
						alleles.addAll(g.getAlleles());
						alleles.remove(Allele.NO_CALL);
						vcb.alleles(alleles);
						vcb.genotypes(Collections.singletonList(gb.make()));
						final VariantContext ctx2= this.recalculator.apply(vcb.make());
						w.add(ctx2);
						}
					}
				for(final String sample:samples)
					{
					LOG.info("Closing for sample "+sample);
					CloserUtil.close(sample2writer.get(sample));
					}
				progress.close();
				out.flush();
				return RETURN_OK;
				}
			catch (final Exception e)
				{
				LOG.error(e);
				return -1;
				}
			finally
				{
				CloserUtil.close(out);
				CloserUtil.close(in);
				}
			}

	public static void main(final String[] args)
		{
		new Biostar130456().instanceMainWithExit(args);
		}
	}
