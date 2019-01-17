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
package com.github.lindenb.jvarkit.tools.vcffixindels;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**

BEGIN_DOC


### See also

 *  https://github.com/lindenb/jvarkit/wiki/VCFFixIndels
 *  "Unified Representation of Genetic Variants" http://bioinformatics.oxfordjournals.org/content/early/2015/02/19/bioinformatics.btv112.abstract (hey ! it was published after I wrote this tool !)
 *  https://github.com/quinlan-lab/vcftidy/blob/master/vcftidy.py
 *  http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/

### Example

```
$ curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/input_callsets/si/ALL.wgs.samtools_pass_filter.20130502.snps_indels.low_coverage.sites.vcf.gz" |\
 gunzip -c | java -jar dist/vcfstripannot.jar -k '*' 2> /dev/null |\
 java -jar dist/vcffixindels.jar  2> /dev/null | grep FIX | head -n 15

##INFO=<ID=INDELFIXED,Number=1,Type=String,Description="Fix Indels for @SolenaLS (position|alleles...)">
1   2030133 .   T   TTTTGT,TTTTG    999 PASS    INDELFIXED=2030101|CGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGT*|CGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGT|CGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTG
1   3046430 .   C   CCCT,CCC    999 PASS    INDELFIXED=3046429|TC*|TCCCT|TCCC
1   4258325 rs137902679;rs61115653  A   AAT,AA  999 PASS    INDELFIXED=4258316|CAAAAAAAAA*|CAAAAAAAAAA|CAAAAAAAAAAT
1   5374885 rs59294415  C   CCCC,CCCCA  999 PASS    INDELFIXED=5374881|TCCCC*|TCCCCCCC|TCCCCCCCA
1   5669438 rs143435517 C   CACAT,CAC   999 PASS    INDELFIXED=5669414|TACACACACACACACACACACACAC*|TACACACACACACACACACACACACAC|TACACACACACACACACACACACACACAT
1   5702062 .   A   AA,AAC  999 PASS    INDELFIXED=5702060|TAA*|TAAAC|TAAA
1   5713682 rs70977965  A   AAAAA,AAAAAC    999 PASS    INDELFIXED=5713678|CAAAA*|CAAAAAAAA|CAAAAAAAAC
1   5911136 .   T   TGCCATT,TGCCATTCCAAAGAGGCACTCA  999 PASS    INDELFIXED=5911135|CT*|CTGCCATTCCAAAGAGGCACTCA|CTGCCATT
1   6067269 rs34064079;rs59468731   G   GG,GGC  999 PASS    INDELFIXED=6067261|TGGGGGGGG*|TGGGGGGGGG|TGGGGGGGGGC
1   6069948 .   TC  T,TTC   999 PASS    INDELFIXED=6069933|CTTTTTTTTTTTTTTTC*|CTTTTTTTTTTTTTTTTC|CTTTTTTTTTTTTTTT
1   6480784 .   C   CGGGCCCCAGGCTGCCCGCC,CGGGCCCCAGGCTGCCCGCCT  999 PASS    INDELFIXED=6480783|GC*|GCGGGCCCCAGGCTGCCCGCCT|GCGGGCCCCAGGCTGCCCGCC
1   6829081 rs34184977;rs5772255    A   AAC,AA  999 PASS    INDELFIXED=6829070|TAAAAAAAAAAA*|TAAAAAAAAAAAA|TAAAAAAAAAAAAC
1   7086193 .   AG  A,AAG   999 PASS    INDELFIXED=7086179|TAAAAAAAAAAAAAAG*|TAAAAAAAAAAAAAAAG|TAAAAAAAAAAAAAA
1   8096161 .   T   TATATATATAC,TAT 999 PASS    INDELFIXED=8096143|CATATATATATATATATAT*|CATATATATATATATATATAT|CATATATATATATATATATATATATATAC

```


END_DOC
*/
@Program(
	name="vcffixindels",
	description="Fix samtools indels (for @SolenaLS)",
	keywords={"vcf","indel"},
	deprecatedMsg="use `bcftools norm`"
	)
public class VCFFixIndels extends Launcher
	{
	private static final Logger LOG = Logger.build(VCFFixIndels.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-t","--tag"},description="INFO/tag")
	private String tag = "INDELFIXED";

	public VCFFixIndels()
		{
		}
	
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator r,
			final VariantContextWriter w
			)
		{
		long nChanged=0L;
		final VCFHeader header=r.getHeader();
		
		final VCFHeader h2=new VCFHeader(header);
		addMetaData(h2);
		h2.addMetaDataLine(new VCFInfoHeaderLine(this.tag,1,VCFHeaderLineType.String,"Fix Indels for @SolenaLS (position|alleles...)"));
		JVarkitVersion.getInstance().addMetaData(this, h2);
		w.writeHeader(h2);
		final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
		
		while(r.hasNext())
			{
			boolean somethingChanged=false;
			final VariantContext ctx=progress.watch(r.next());
			/* map old allele to new allele */
			final Map<Allele, Allele> original2modified=new HashMap<Allele, Allele>();
			/* did we found a strange allele (symbolic, etc ) */
			boolean strange_allele_flag=false;
			for(final Allele a:ctx.getAlleles())
				{
				original2modified.put(a, a);
				if(a.isSymbolic() || a.isNoCall() || a.equals(Allele.SPAN_DEL))
					{
					strange_allele_flag=true;
					break;
					}
				}
			
			if(strange_allele_flag ||
				original2modified.size()<2 /* at least 2 alleles: REF+ ALT */
				)
				{
				w.add(ctx);
				continue;
				}
			
			/* record chromStart if we need to shift to the right */
			int chromStart= ctx.getStart();
			
			/* trim right then left  */
			for(int side=0;side<2;++side)
				{
				boolean done=false;
				while(!done)
					{
					boolean ok_side=true;
					/* common nucleotide at end/start */
					Character targetChar = null;
					done=true;
					//scan side
					Set<Allele> keys = new HashSet<>(original2modified.keySet());  
					for(Allele k:keys)
						{
						Allele newAllele = original2modified.get(k);
						if(newAllele.isSymbolic()) 
							{
							ok_side = false;
							break;
							}
						String baseString = newAllele.getBaseString().trim().toUpperCase();
						if(baseString.length()<2)
							{
							ok_side = false;
							break;
							}
						/* first or last char or all sequences
						 * side==0 : right
						 * side==1 : left
						 * */
						Character baseChar=
							(side==0?
							baseString.charAt(baseString.length()-1):
							baseString.charAt(0)
							);
						if(targetChar==null)
							{
							targetChar = baseChar;
							}
						else if(!targetChar.equals(baseChar))
							{
							/* doesn't end with same nucleotide */
							ok_side = false;
							break;
							}
						}
					/* ok we can shift all alleles */
					if(ok_side && targetChar!=null)
						{
						done=false;
						somethingChanged=true;
						for(final Allele k:keys)
							{
							Allele newAllele = original2modified.get(k);
							String baseString = newAllele.getBaseString().trim().toUpperCase();
							if(side==0)//remove last nucleotide
								{
								newAllele = Allele.create(
									baseString.substring(0,baseString.length()-1),
									newAllele.isReference());
								}
							else
								{
								newAllele = Allele.create(
										baseString.substring(1),
										newAllele.isReference());
								
								}
							original2modified.put(k, newAllele);
							}
						if(side==1) chromStart++;
						}
					}/* end of while done */
				
				}/* end side */
			
			if(!somethingChanged)
				{
				w.add(ctx);
				continue;
				}
			
			final VariantContextBuilder b=new VariantContextBuilder(ctx);
			b.start(chromStart);
			final Allele newRef=original2modified.get(ctx.getReference());
			b.stop(chromStart+newRef.getBaseString().length()-1);
			b.alleles(original2modified.values());
			final List<Genotype> genotypes=new ArrayList<>();
			for(final String sample:header.getSampleNamesInOrder())
				{
				final Genotype g = ctx.getGenotype(sample);
				if(!g.isCalled()) 
					{
					genotypes.add(g);
					continue;
					}
				final GenotypeBuilder gb=new GenotypeBuilder(g);
				final List<Allele> aL=new ArrayList<>();
				for(final Allele a:g.getAlleles())
					{
					aL.add(original2modified.get(a));
					}
				gb.alleles(aL);
				genotypes.add(gb.make());
				}
			
			final StringBuilder tagContent=new StringBuilder();
			tagContent.append(String.valueOf(ctx.getStart()));
			for(final Allele a:ctx.getAlleles())
				{
				tagContent.append("|");
				tagContent.append(a.toString());
				}
			b.attribute(this.tag,tagContent.toString());
			b.genotypes(genotypes);
			
			w.add( b.make());
			++nChanged;
			if(w.checkError()) break;
			}
		progress.finish();
		LOG.info("indels changed:"+nChanged);
		return RETURN_OK;
		}

	
	@Override
	public int doWork(final List<String> args) {
		if(StringUtil.isBlank(this.tag)) {
			LOG.error("empty INFO/tag");
			return -1;
		}
		try {
			return doVcfToVcf(args,outputFile);
			} 
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		}
	
	
	public static void main(final String[] args) throws IOException
		{
		new VCFFixIndels().instanceMainWithExit(args);
		}
}
