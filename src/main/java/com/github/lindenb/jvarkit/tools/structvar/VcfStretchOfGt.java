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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/vcfstrechofgt.jar -p src/test/resources/test_vcf01.ped src/test/resources/test_vcf01.vcf

#chrom	start0	end0	length	count.affected.variants	average.affected.depth	count.other.variants
1	870316	870317	1	1	5.0	0
1	919500	919501	1	1	0.0	0
1	963703	963704	1	1	0.0	0
1	1004201	1004202	1	1	0.0	0
```

END_DOC
*/
@Program(name="vcfstrechofgt",
description="Try to finds deletion by searching strech of HOM_REF/HOM_VAR/NO_CALL Genotypes.",
keywords={"vcf","deletion"}
)
public class VcfStretchOfGt extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfStretchOfGt.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-p","--pedigree"},description=Pedigree.OPT_DESCRIPTION)
	private File pedigreeFile = null;
	
	private class Stretch
		{
		String contig;
		int start;
		int end;
		int countVariants=0;
		double sumAvgDp=0.0;
		int countOthers=0;
		}
	private void dump(final PrintWriter w,final Stretch data) {
		if(data.countOthers>=data.countVariants) return;
		w.print(data.contig);
		w.print('\t');
		w.print(data.start-1);
		w.print('\t');
		w.print(data.end);
		w.print('\t');
		w.print(1+data.end-data.start);
		w.print('\t');
		w.print(data.countVariants);
		w.print('\t');
		w.print(data.sumAvgDp/data.countVariants);
		w.print('\t');
		w.print(data.countOthers);
		w.println();
	}
	
	@Override
	public int doWork(final List<String> args) {
		VCFIterator iter=null;
		PrintWriter w=null;
		try {
			iter = super.openVCFIterator(oneFileOrNull(args));
			final VCFHeader header=iter.getHeader();
			if(!header.hasGenotypingData()) {
				LOG.error("No genotype in input");
				return -1;
			}
			final Pedigree pedigree = (this.pedigreeFile==null?
					Pedigree.newParser().parse(header):
					Pedigree.newParser().parse(this.pedigreeFile)
					);
			pedigree.verifyPersonsHaveUniqueNames();
			final Set<String> affected = pedigree.getAffected().
					stream().
					map(P->P.getId()).
					filter(ID->header.getSampleNameToOffset().containsKey(ID)).
					collect(Collectors.toSet());
			final Set<String> otherSamples = new HashSet<>(header.getSampleNamesInOrder());
			otherSamples.removeAll(affected);
			
			if(affected.isEmpty()) {
				LOG.error("No affected sample in VCF");
				return -1;
				}
			w = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			w.print("#chrom");
			w.print('\t');
			w.print("start0");
			w.print('\t');
			w.print("end0");
			w.print('\t');
			w.print("length");
			w.print('\t');
			w.print("count.affected.variants");
			w.print('\t');
			w.print("average.affected.depth");
			w.print('\t');
			w.print("count.other.variants");
			w.println();
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			Stretch current=null;
			
			final Predicate<Genotype> acceptGenotype=G->G.isHomRef() || G.isHomVar()|| G.isNoCall();
			
			while(iter.hasNext()) {
				final VariantContext ctx=progress.apply(iter.next());
				if(current!=null && !current.contig.equals(ctx.getContig())) {
					dump(w,current);
					current=null;
					}
				boolean ok=
						affected.stream().
						map(S->ctx.getGenotype(S)).
						allMatch(acceptGenotype);
				
				boolean otherOk=
						otherSamples.stream().
						map(S->ctx.getGenotype(S)).
						allMatch(acceptGenotype);
				
				if(!ok) {
					if(current!=null) dump(w,current);
					current=null;
					}
				else if(current==null) {
					current=new Stretch();
					current.contig=ctx.getContig();
					current.start= ctx.getStart();
					}
				
				//do this for any Stretch
				if(current!=null) {
					current.end= ctx.getEnd();
					current.countVariants++;
					current.sumAvgDp += ctx.getGenotypes().stream().
							filter(G->G.hasDP() && affected.contains(G.getSampleName())).
							mapToInt(G->G.getDP()).
							average().
							orElse(0.0);
					current.countOthers += (otherOk?1:0);
					}
				}
			if(current!=null) dump(w,current);
			w.flush();w.close();w=null;
			iter.close();iter=null;
			progress.close();
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(iter);
			CloserUtil.close(w);
			}
		}

	public static void main(final String[] args)
		{
		new VcfStretchOfGt().instanceMainWithExit(args);
		}	

}
