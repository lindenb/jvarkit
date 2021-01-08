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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;

import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
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

## See also

```
bcftools roh
```

END_DOC
*/
@Program(name="vcfstrechofgt",
description="Try to finds deletion by searching strech of HOM_REF/HOM_VAR/NO_CALL Genotypes.",
keywords={"vcf","deletion","cnv"},
creationDate="20190103",
modificationDate="20200108"
)
public class VcfStretchOfGt extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfStretchOfGt.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--pedigree"},description="If defined, the tool will use the affected sample and find strech where all affected could be a DEL." + PedigreeParser.OPT_DESC)
	private Path pedigreePath = null;
	@Parameter(names={"-a","--affected"},description="Same as option --pedigree but provide the name of the samples using a comma/space/semicolon separated string")
	private String affectedStr = null;
	@Parameter(names={"-nc","--no-call"},description="Do not accept NO_CALL genotypes.")
	private boolean exclude_no_call = false;
	
	
	
	private static class Stretch
		{
		String contig;
		int start;
		int end;
		int countVariants=0;
		double sumAvgDp=0.0;
		int countOthers=0;
		}

	
	private class SampleSet {
		final Set<String> affected;
		final Set<String> otherSamples;
		Stretch current = null;
		final Predicate<Genotype> acceptGenotype;
		
		SampleSet(final Set<String> all_samples,final Set<String> affected) {
			this.affected = affected;
			this.otherSamples = new HashSet<>(all_samples);
			this.otherSamples.removeAll(this.affected);
			this.acceptGenotype=G->G.isHomRef() || G.isHomVar()|| (G.isNoCall() && !exclude_no_call);
			}
		
		private void dump(final PrintWriter w) {
			if(current==null) return ;
			if(current.countOthers>=current.countVariants) return;
			w.print(current.contig);
			w.print('\t');
			w.print(current.start-1);
			w.print('\t');
			w.print(current.end);
			w.print('\t');
			w.print(CoordMath.getLength(current.start, current.end));
			w.print('\t');
			w.print(String.join(";",this.affected));
			w.print('\t');
			w.print(current.countVariants);
			w.print('\t');
			w.print(String.format("%.2f",current.sumAvgDp/current.countVariants));
			w.print('\t');
			w.print(current.countOthers);
			w.println();

		}
		
		private void visit(final PrintWriter w,final VariantContext ctx) {
			if(current!=null && !current.contig.equals(ctx.getContig())) {
				dump(w);
				current=null;
				}
			final boolean ok=
					affected.stream().
					map(S->ctx.getGenotype(S)).
					allMatch(acceptGenotype);
			
			final boolean otherOk=
					otherSamples.stream().
					map(S->ctx.getGenotype(S)).
					allMatch(acceptGenotype);
			
			
			
			if(!ok) {
				if(current!=null) dump(w);
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
			final Set<String> vcf_samples = new HashSet<>(header.getSampleNamesInOrder());
			
			final List<SampleSet> all_sample_set = new ArrayList<>();
			
			if(!StringUtils.isBlank(this.affectedStr)) {
				
				final Set<String> affected = Arrays.stream(this.affectedStr.split("[ ,;]+")).
						filter(S->!StringUtils.isBlank(S)).
						filter(S->vcf_samples.contains(S)).
						collect(Collectors.toSet());
				
				if(affected.isEmpty()) {
					LOG.error("No affected sample in string "+this.affectedStr);
					return -1;
					}
				
				/* at least 2 because singleton sample are created later */
				if(affected.size()>1) {
					all_sample_set.add(new SampleSet(vcf_samples,affected));
					}
				}
			
			if(this.pedigreePath!=null) {
				final Pedigree pedigree =  new PedigreeParser().parse(this.pedigreePath);
				
				final Set<String> affected = pedigree.
						getAffectedSamples().
						stream().
						map(P->P.getId()).
						filter(ID->header.getSampleNameToOffset().containsKey(ID)).
						collect(Collectors.toSet());
				
				if(affected.isEmpty()) {
					LOG.error("No affected sample in pedigree "+this.pedigreePath);
					return -1;
					}
				
				/* at least 2 because singleton sample are created later */
				if(affected.size()>1) {
					all_sample_set.add(new SampleSet(vcf_samples,affected));
					}
				}
			
			/* singletons */
			for(final String sn:vcf_samples) {
				final SampleSet sampleSet = new SampleSet(vcf_samples,Collections.singleton(sn));
				all_sample_set.add(sampleSet);
			}
			
		
			
			
			w = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			
			w.print("#chrom");
			w.print('\t');
			w.print("start0");
			w.print('\t');
			w.print("end0");
			w.print('\t');
			w.print("length");
			w.print('\t');
			w.print("samples");
			w.print('\t');
			w.print("count.affected.variants");
			w.print('\t');
			w.print("average.affected.depth");
			w.print('\t');
			w.print("count.other.variants");
			w.println();
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			
			
			while(iter.hasNext()) {
				final VariantContext ctx=progress.apply(iter.next());
				for(final SampleSet snSet:all_sample_set) snSet.visit(w, ctx);
				}
			for(final SampleSet snSet:all_sample_set) snSet.dump(w);
			w.flush();w.close();w=null;
			iter.close();iter=null;
			progress.close();
			return 0;
			}
		catch(final Throwable err) {
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
