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
/**

BEGIN_DOC




```

 curl -s  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/integrated_call_sets/ALL.chr22.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz" | gunzip -c | \
java -jar dist/allelefreqcalc.jar | head

CHR	POS	ID	REF	ALT	TOTAL_CNT	ALT_CNT	FRQ
22	16050408	rs149201999	T	C	2184	134	0.06135531
22	16050612	rs146752890	C	G	2184	184	0.08424909
22	16050678	rs139377059	C	T	2184	113	0.051739927
22	16050984	rs188945759	C	G	2184	5	0.0022893774
22	16051107	rs6518357	C	A	2184	127	0.058150183
22	16051249	rs62224609	T	C	2184	157	0.07188645
22	16051347	rs62224610	G	C	2184	650	0.29761904
22	16051453	rs143503259	A	C	2184	160	0.07326008
22	16051477	rs192339082	C	A	2184	2	9.157509E-4

```




END_DOC
*/


import java.io.File;
import java.io.PrintStream;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

@Program(
	name="allelefreqcalc",
	description="Allele Frequency Calculator",
	keywords={"vcf","af"},
	deprecatedMsg="Use bioalcidae"
	)
public class AlleleFrequencyCalculator extends Launcher
	{
	private static final Logger LOG = Logger.build(AlleleFrequencyCalculator.class).make();
	@Parameter(names={"-o","--output"},description=Launcher.OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	public AlleleFrequencyCalculator()
		{
		}
	@Override
	public int doWork(final List<String> args)
		{
		PrintStream out = null;
		VCFIterator in = null;
		try
			{
			in = super.openVCFIterator(oneAndOnlyOneFile(args));
		
			out = openFileOrStdoutAsPrintStream(outputFile);
			
			
			out.println("CHR\tPOS\tID\tREF\tALT\tTOTAL_CNT\tALT_CNT\tFRQ");
			while(in.hasNext() && !out.checkError())
				{
				
				final VariantContext ctx=in.next();
				final Allele ref=ctx.getReference();
				if(ref==null) continue;
				if(ctx.getNSamples()==0 || ctx.getAlternateAlleles().isEmpty()) continue;
				final Allele alt=ctx.getAltAlleleWithHighestAlleleCount();
				if(alt==null) continue;
				
				final GenotypesContext genotypes=ctx.getGenotypes();
				if(genotypes==null) continue;
				int total_ctn=0;
				int alt_ctn=0;
				for(int i=0;i< genotypes.size();++i)
					{
					final Genotype g=genotypes.get(i);
					for(final Allele allele: g.getAlleles())
						{
						if(allele.equals(ref))
							{
							total_ctn++;
							}
						else if (allele.equals(alt))
							{
							total_ctn++;
							alt_ctn++;
							}
						}
					
					}
				
				
				out.print(ctx.getContig());
				out.print("\t");
				out.print(ctx.getStart());
				out.print("\t");
				out.print(ctx.hasID()?ctx.getID():".");
				out.print("\t");
				out.print(ref.getBaseString());
				out.print("\t");
				out.print(alt.getBaseString());
				out.print("\t");
				out.print(total_ctn);
				out.print("\t");
				out.print(alt_ctn);
				out.print("\t");
				out.print(alt_ctn/(float)total_ctn);
				out.println();
				}
			out.flush();
			out.close();
			out = null;
			
		return RETURN_OK;
		}
	catch(final Exception err)
		{
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(out);
		CloserUtil.close(in);
		}
	}
	 	

			 
	
	public static void main(String[] args) {
		new AlleleFrequencyCalculator().instanceMainWithExit(args);

	}

}
