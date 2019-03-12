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
package com.github.lindenb.jvarkit.tools.vcfpedcheckrare;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Trio;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory.AFExtractor;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/*
BEGIN_DOC

 
## Example


```bas
```

END_DOC
 */
@Program(name="vcfpedcheckrare",
	description="Check pedigree transmission by looking at the rare variants",
	deprecatedMsg="use GATK combineVariants.",
	keywords={"vcf","pedigree","rare","sample","trio"}
	)
public class VcfPedigreeCheckRare extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfPedigreeCheckRare.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--ped","--pedigree"},description=PedigreeParser.OPT_DESC,required=true)
	private Path pedFile = null;

	private enum Type {
		INHERITED,DENOVO
		}
	
	@Override
	public int doWork(final List<String> args) {
		VCFIterator iter=null;
		try {
			iter = super.openVCFIterator(oneFileOrNull(args));
			final VCFHeader header = iter.getHeader();
			if(!header.hasGenotypingData()) {
				LOG.error("no genotype data in input.");
				return -1;	
				}
			
			
			final Set<String> sampleSet = new HashSet<>(header.getSampleNamesInOrder());
			
			final Pedigree pedigree = new PedigreeParser().parse(this.pedFile);
			final TreeSet<Trio> trios = pedigree.getTrios().
				stream().
				filter(T->sampleSet.contains(T.getChild().getId())).
				filter(T->(T.hasFather() && sampleSet.contains(T.getFather().getId())) ||
						  (T.hasMother() && sampleSet.contains(T.getMother().getId()))
						  ).
				collect(Collectors.toCollection(TreeSet::new));
			
			if(trios.isEmpty())
				{
				LOG.info("no trio found with a correspondance in the VCF header");
				return -1;
				}
			final Map<String,Counter<Type>> sample2count = new TreeMap<>();

			for(final Trio trio:trios) {
				sample2count.put(trio.getChild().getId(), new Counter<>());
				}
				
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			while(iter.hasNext()) 
				{
				final VariantContext ctx  = progress.apply(iter.next());
				if(!ctx.isVariant()) continue;
				final List<Allele> alt_alleles = ctx.getAlternateAlleles();
				
				
				
				for(int alleleidx=0 ; alleleidx < alt_alleles.size(); ++alleleidx) {
					final Allele alt = alt_alleles.get(alleleidx);
					
					/** loop over trios */
					for(final Trio trio: trios) {
						final Genotype gc = ctx.getGenotype(trio.getChild().getId());
						/* child must be HET (RARE) */
						if(gc==null || !gc.isHet()) continue;
						/* child must carry alt allele */
						if(!gc.getAlleles().contains(alt)) continue;
						
						
						final Genotype gf = trio.hasFather()?null:ctx.getGenotype(trio.getFather().getId());
						final Genotype gm = trio.hasMother()?null:ctx.getGenotype(trio.getMother().getId());
						
						final boolean parent_have_alt=
								(gf!=null && gf.isHet() && gf.getAlleles().contains(alt)) ||
								(gm!=null && gm.isHet() && gm.getAlleles().contains(alt))
								;
						final Type type=parent_have_alt?
								Type.INHERITED:
								Type.DENOVO
								;
						sample2count.get(gc.getSampleName()).incr(type);
						}
					}
				
				}
			
			progress.close();
			
			/* done, write table */
			final PrintWriter pw = super.openPathOrStdoutAsPrintWriter(outputFile);
			pw.println("#sample");
			for(final Type f : Type.values()) pw.print("\t"+f.name());
			pw.println();
			for(final String sample: sample2count.keySet())
				{
				final Counter<Type> count = sample2count.get(sample);
				pw.print(sample);
				for(final Type f : Type.values()) {
					pw.print("\t");
					pw.print(count.count(f));
					}
				pw.println();
				}
				
			pw.flush();
			pw.close();
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	
	public static void main(final String[] args) {
		new VcfPedigreeCheckRare().instanceMainWithExit(args);
		}
	
}
