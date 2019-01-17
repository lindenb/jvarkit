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
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.script.Bindings;
import javax.script.CompiledScript;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;

/**

BEGIN_DOC




### Example

The script injects in the context:
 *  header a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html
 *  variant a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html
 *  genotype a https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Genotype.html

if the returned value is false, the genotype is set to no call or to hom-ref.


```
$ cat  ~/src/gatk-ui/testdata/mutations.vcf | java -jar dist/vcfresetgenotypejs.jar -homref -e '!genotype.isHomRef()' |grep -v "##"

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  S1      S2      S3      S4
rotavirus       51      .       A       G       22.55   .       AC1=2;AF1=0.25;BQB=1;DP=944;DP4=849,0,93,0;FQ=23.7972;G3=0.75,0,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.993129;SGB=-61.9012;VDB=3.53678e-05  GT:PL   0/0:0,255,134   0/0:0,255,127   0/0:0,255,137   1/1:70,255,0
rotavirus       91      .       A       T       5.45    .       AC1=1;AF1=0.124963;BQB=0.951201;DP=1359;DP4=1134,0,225,0;FQ=5.8713;MQ=60;MQ0F=0;MQB=1;PV4=1,4.80825e-05,1,1;RPB=0.0393173;SGB=-369.163;VDB=0.313337 GT:PL   0/0:0,255,133   0/1:40,0,31     0/0:0,255,134   0/0:0,255,82
rotavirus       130     .       T       C       4.12    .       AC1=1;AF1=0.124933;BQB=1;DP=1349;DP4=1139,0,204,0;FQ=4.48321;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.762964;SGB=-335.275;VDB=0.00084636  GT:PL0/1:38,0,35      0/0:0,255,132   0/0:0,255,132   0/0:0,255,79


```


 */
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
 BEGIN_DOC
 
 The script injects in the context:

  * **header** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/vcf/VCFHeader.html)
  * **variant** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html)
  * **genotype** a [https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Genotype.html](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/Genotype.html)



if the returned value is false, the genotype is set to no call or to hom-ref.

```bash
$ cat  ~/src/gatk-ui/testdata/mutations.vcf | java -jar dist/vcfresetgenotypejs.jar -homref -e '!genotype.isHomRef()' |grep -v "##"

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC1=2;AF1=0.25;BQB=1;DP=944;DP4=849,0,93,0;FQ=23.7972;G3=0.75,0,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.993129;SGB=-61.9012;VDB=3.53678e-05	GT:PL	0/0:0,255,134	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	AC1=1;AF1=0.124963;BQB=0.951201;DP=1359;DP4=1134,0,225,0;FQ=5.8713;MQ=60;MQ0F=0;MQB=1;PV4=1,4.80825e-05,1,1;RPB=0.0393173;SGB=-369.163;VDB=0.313337	GT:PL	0/0:0,255,133	0/1:40,0,31	0/0:0,255,134	0/0:0,255,82
rotavirus	130	.	T	C	4.12	.	AC1=1;AF1=0.124933;BQB=1;DP=1349;DP4=1139,0,204,0;FQ=4.48321;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.762964;SGB=-335.275;VDB=0.00084636	GT:PL	0/1:38,0,35	0/0:0,255,132	0/0:0,255,132	0/0:0,255,79
```



 
 END_DOC
 */
@Program(name="vcfremovegenotypejs",description="Reset Genotype in VCF using a javascript expression")
public class VcfRemoveGenotypeJs extends Launcher {
	private static final Logger LOG = Logger.build(VcfRemoveGenotypeJs.class).make();


	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;
	
	@Parameter(names={"-e","--expression"},description=" (js expression). Optional.")
	private String scriptExpr=null;
	@Parameter(names={"-f","--script"},description=" (js file). Optional.")
	private File scriptFile=null;


	@Parameter(names={"-g","--removeCtxNoGenotype"},description="Remove variants having no called genotype or all home Ref.")
	private boolean removeCtxNoGenotype = false;

	@Parameter(names={"-homref","--homref"},description="Replace variant with homref instead of nocall")
	private boolean replaceByHomRef = false;

	@Parameter(names={"-filter","--filter"},description="if not empty, don't delete the genotype but filter it.")
	private String filterName = "";
	
	private CompiledScript  script=null;

	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter out) {
		try {
		this.script  = super.compileJavascript(scriptExpr, scriptFile);
		final	VCFHeader h2=new VCFHeader(in.getHeader());
		
		if(!this.filterName.isEmpty() )
			{
			h2.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY, 1, VCFHeaderLineType.String, "Genotype-level filter"));
			}	
		addMetaData(h2);
		out.writeHeader(h2);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(in.getHeader());
        final Bindings bindings = this.script.getEngine().createBindings();
        bindings.put("header", in.getHeader());

		while(in.hasNext())
			{
			final VariantContext ctx = progress.watch(in.next());
			bindings.put("variant", ctx);
			final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
			List<Genotype> genotypes= new ArrayList<>();
			int countCalled = ctx.getNSamples();
			
			for(int i=0;i< ctx.getNSamples();++i)
				{
				Genotype genotype = ctx.getGenotype(i);
				bindings.put("genotype", genotype);
				
				if(!genotype.isCalled() || genotype.isNoCall() || !genotype.isAvailable())
					{
					countCalled--;
					}
				else if(genotype.isCalled() &&
					!super.evalJavaScriptBoolean(this.script, bindings)) {
					if(!this.filterName.isEmpty()) {
						if(!genotype.isFiltered()) 
							{
							genotype = new GenotypeBuilder(genotype).filters(this.filterName).make();
							}
						}
					else if(this.replaceByHomRef){
						List<Allele> homRefList=new ArrayList<>(genotype.getPloidy());
						for(int p=0;p< genotype.getPloidy();++p)
							{
							homRefList.add(ctx.getReference());
							}
						genotype = new GenotypeBuilder(genotype).alleles(homRefList).make();
						} 
					else
						{
						genotype = GenotypeBuilder.createMissing(genotype.getSampleName(), genotype.getPloidy());
						}
					countCalled--;
					}
				genotypes.add(genotype);
				}
			if(countCalled==0 && this.removeCtxNoGenotype) {
				continue;
			}
			vcb.genotypes(genotypes);
			out.add(vcb.make());
			}
		
		progress.finish();
		
		return RETURN_OK;
		} catch(Exception err) {
			LOG.error(err);
			return -1;
		} finally
		{
			this.script=null;
		}
		}
	
	@Override
	public int doWork(List<String> args) {
		return doVcfToVcf(args, outputFile);
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new VcfRemoveGenotypeJs().instanceMainWithExit(args);
	}

}
