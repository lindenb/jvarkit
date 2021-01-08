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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

/**
BEGIN_DOC

## Example

```
t$ find /home/lindenb/src/jvarkit/src/test/resources -name "*.vcf.gz" -o -name "*.bam" -o -name "*.cram" | java -jar dist/howmanybamdict.jar | cut -c 1-${COLUMNS}
9a5c58c2c91e731135b27ed14974523a	.	85	3101976562	1=249250621;2=243199373;3=198022430;4=191154276;5=180915260;6=17111
9a5c58c2c91e731135b27ed14974523a	/home/lindenb/src/jvarkit/src/test/resources/ExAC.r1.sites.vep.vcf.gz
4677ece43eea2b029d0d33fe130ea6c7	.	86	3137454505	chr1=249250621;chr2=243199373;chr3=198022430;chr4=191154276;chr5=18
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/manta.B00I9CJ.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	.	11	18490	RF01=3302;RF02=2687;RF03=2592;RF04=2362;RF05=1579;RF06=1356;RF07=1074;RF
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S2.vcf.gz
9a5c58c2c91e731135b27ed14974523a	/home/lindenb/src/jvarkit/src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S3.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S1.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/rotavirus_rf.unifiedgenotyper.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S4.vcf.gz
635de5cb51973d45844fa713ac0b7719	.	2	85	ref=45;ref2=40	/home/lindenb/src/jvarkit/src/test/resources/toy.vcf.gz
635de5cb51973d45844fa713ac0b7719	/home/lindenb/src/jvarkit/src/test/resources/toy.vcf.gz
635de5cb51973d45844fa713ac0b7719	/home/lindenb/src/jvarkit/src/test/resources/toy.bam
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S4.bam
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/roxan.hs37d5.csq.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/rotavirus_rf.ann.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/rotavirus_rf.vcf.gz
df8200dd2a49e25bc98df5f2c45ac36a	.	2	85	ref=45;ref2=40	/home/lindenb/src/jvarkit/src/test/resources/toy.cram
df8200dd2a49e25bc98df5f2c45ac36a	/home/lindenb/src/jvarkit/src/test/resources/toy.cram
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/rotavirus_rf.freebayes.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S1.bam
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S2.bam
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S5.bam
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/manta.B00GWIU.vcf.gz
5de1aef8e9e3be34ae7696779a4791c7	.	3366	3217346917	chr1=248956422;chr2=242193529;chr3=198295559;chr4=190214555;chr5=
5de1aef8e9e3be34ae7696779a4791c7	/home/lindenb/src/jvarkit/src/test/resources/FAB23716.nanopore.bam
e07f3b093938833945fa357c4b37bdf9	.	292	3100014256	chr1=248956422;chr2=242193529;chr3=198295559;chr4=190214555;chr5=1
e07f3b093938833945fa357c4b37bdf9	/home/lindenb/src/jvarkit/src/test/resources/ENCFF331CGL.rnaseq.b38.bam
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/retrocopy01.bwa.bam
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/manta.B00GWGD.vcf.gz
4677ece43eea2b029d0d33fe130ea6c7	/home/lindenb/src/jvarkit/src/test/resources/manta.D000Q1R.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S5.vcf.gz
bd7e0928fc3c810e48fafc53a4222ed5	/home/lindenb/src/jvarkit/src/test/resources/S3.bam
f8d942cb3fc6ebef618a0b0ba3f4ef99	.	24	3095677412	1=249250621;2=243199373;3=198022430;4=191154276;5=180915260;6=17111
f8d942cb3fc6ebef618a0b0ba3f4ef99	/home/lindenb/src/jvarkit/src/test/resources/gnomad_v2_sv.sites.vcf.gz
4327878bf3a3073edf6e77fda48033fe	.	86	3137454505	1=249250621;2=243199373;3=198022430;4=191154276;5=180915260;6=17111
4327878bf3a3073edf6e77fda48033fe	/home/lindenb/src/jvarkit/src/test/resources/HG02260.transloc.chr9.14.bam
9a5c58c2c91e731135b27ed14974523a	/home/lindenb/src/jvarkit/src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz
(...)
```

END_DOC
 */
@Program(name="howmanybamdict",
	description="finds if there's are some differences in the sequence dictionaries.",
	keywords={"sam","bam","dict","vcf"},
	biostars=468541,
	creationDate="20131108",
	modificationDate="20201021"
	)
public class HowManyBamDict extends Launcher {
	private static final Logger LOG = Logger.build(HowManyBamDict.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;

	private final Set<String>  all_md5 = new HashSet<>();
	
 	private void handle(final PrintWriter out,final Path f)
 		{
		final SAMSequenceDictionary dict0;
 		try {
			dict0 = SAMSequenceDictionaryExtractor.extractDictionary(f);
 			}
		catch(final Throwable err) {
			LOG.warn(err);
			out.print("#ERROR");
			out.print("\t");
			out.print(f);
			out.println();
			return;
			}
		
		if(dict0==null) {
			out.print("#DICT_MISSING");
			out.print("\t");
			out.print(f);
			out.println();
			}
		else
			{
			final String md5 = dict0.md5();
			if(!this.all_md5.contains(md5)) {
				this.all_md5.add(md5);
				out.print(md5);
				out.print("\t");
				out.print(".");
				out.print("\t");
				out.print(dict0.size());
				out.print("\t");
				out.print(dict0.getReferenceLength());
				out.print("\t");
				out.print(dict0.getSequences().stream().map(ssr->ssr.getSequenceName()+"="+ssr.getSequenceLength()).collect(Collectors.joining(";")));
				out.print("\t");
				out.print(f);
				out.println();
				}
			out.print(md5);
			out.print("\t");
			out.print(f);
			out.println();
			}
 		}
	
 	@Override
 	public int doWork(List<String> args) {
		PrintWriter out=null;
		try
			{
			out = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			if(args.isEmpty())
				{
				final PrintWriter fout = out;
				try(BufferedReader in=new BufferedReader(new InputStreamReader(stdin()))) {
					in.lines().
						filter(line->!(StringUtils.isBlank(line) || line.endsWith(File.separator) || line.startsWith("#"))).
						forEach(line->handle(fout,Paths.get(line)));
					}
				}
			else
				{
				for(String filename:args)
					{
					handle(out,Paths.get(filename));
					}
				}
			out.flush();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(out);
			}
		}

 	public static void main(final String[] args)
		{
		new HowManyBamDict().instanceMainWithExit(args);
		}

}
