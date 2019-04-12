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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.variant.vcf.VCFIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC

## Deprecated

use `bcftools annotate` with `--rename-chrs file`

> rename chromosomes according to the map in file, with "old_name new_name\n" pairs separated by whitespaces, each on a separate line. 

## Example

```bash
$ curl  "https://raw.github.com/arq5x/gemini/master/test/test5.vep.snpeff.vcf" |\
  java -jar dist/vcfrenamechr.jar -i -C -f src/main/resources/chromnames/hg19_to_g1kv37.tsv |\
cut -f 1-5

(...)
##contig=<ID=1,length=249250621>
##contig=<ID=10,length=135534747>
##contig=<ID=11,length=135006516>
##contig=<ID=GL000202.1,length=40103>
##contig=<ID=12,length=133851895>
##contig=<ID=13,length=115169878>
##contig=<ID=14,length=107349540>
##contig=<ID=15,length=102531392>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=GL000203.1,length=37498>
##contig=<ID=GL000204.1,length=81310>
##contig=<ID=GL000205.1,length=174588>
##contig=<ID=GL000206.1,length=41001>
##contig=<ID=18,length=78077248>
##contig=<ID=GL000207.1,length=4262>
##contig=<ID=19,length=59128983>
##contig=<ID=GL000208.1,length=92689>
##contig=<ID=GL000209.1,length=159169>
##contig=<ID=GL000191.1,length=106433>
##contig=<ID=GL000192.1,length=547496>
##contig=<ID=2,length=243199373>
##contig=<ID=20,length=63025520>
##contig=<ID=21,length=48129895>
##contig=<ID=GL000210.1,length=27682>
##contig=<ID=22,length=51304566>
##contig=<ID=3,length=198022430>
##contig=<ID=4,length=191154276>
##contig=<ID=GL000193.1,length=189789>
##contig=<ID=GL000194.1,length=191469>
##contig=<ID=5,length=180915260>
##contig=<ID=6,length=171115067>
##contig=<ID=7,length=159138663>
##contig=<ID=GL000195.1,length=182896>
##contig=<ID=8,length=146364022>
##contig=<ID=GL000196.1,length=38914>
##contig=<ID=GL000197.1,length=37175>
##contig=<ID=9,length=141213431>
##contig=<ID=GL000198.1,length=90085>
##contig=<ID=GL000199.1,length=169874>
##contig=<ID=GL000200.1,length=187035>
##contig=<ID=GL000201.1,length=36148>
##contig=<ID=GL000211.1,length=166566>
##contig=<ID=GL000212.1,length=186858>
##contig=<ID=GL000213.1,length=164239>
##contig=<ID=GL000214.1,length=137718>
##contig=<ID=GL000215.1,length=172545>
##contig=<ID=GL000216.1,length=172294>
##contig=<ID=GL000217.1,length=172149>
##contig=<ID=GL000218.1,length=161147>
##contig=<ID=GL000219.1,length=179198>
##contig=<ID=GL000220.1,length=161802>
##contig=<ID=GL000221.1,length=155397>
##contig=<ID=GL000222.1,length=186861>
##contig=<ID=GL000223.1,length=180455>
##contig=<ID=GL000224.1,length=179693>
##contig=<ID=GL000225.1,length=211173>
##contig=<ID=GL000226.1,length=15008>
##contig=<ID=GL000227.1,length=128374>
##contig=<ID=GL000228.1,length=129120>
##contig=<ID=GL000229.1,length=19913>
##contig=<ID=GL000230.1,length=43691>
##contig=<ID=GL000231.1,length=27386>
##contig=<ID=GL000232.1,length=40652>
##contig=<ID=GL000233.1,length=45941>
##contig=<ID=GL000234.1,length=40531>
##contig=<ID=GL000235.1,length=34474>
##contig=<ID=GL000236.1,length=41934>
##contig=<ID=GL000237.1,length=45867>
##contig=<ID=GL000238.1,length=39939>
##contig=<ID=GL000239.1,length=33824>
##contig=<ID=GL000240.1,length=41933>
##contig=<ID=GL000241.1,length=42152>
##contig=<ID=GL000242.1,length=43523>
##contig=<ID=GL000243.1,length=43341>
##contig=<ID=GL000244.1,length=39929>
##contig=<ID=GL000245.1,length=36651>
##contig=<ID=GL000246.1,length=38154>
##contig=<ID=GL000247.1,length=36422>
##contig=<ID=GL000248.1,length=39786>
##contig=<ID=GL000249.1,length=38502>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##reference=file:///m/cphg-quinlan/cphg-quinlan/shared/genomes/hg19/bwa/gatk/hg19_gatk.fa
#CHROM	POS	ID	REF	ALT
1	145273345	.	T	C
1	156011444	.	T	C
5	64982321	.	T	C
10	1142208	.	T	C
10	126678092	.	G	A
10	135210791	.	T	C
13	48873835	.	G	A
20	36779424	.	G	A
X	17819377	.	T	C


```

## See also

* https://github.com/lindenb/jvarkit/blob/master/src/main/resources/chromnames/g1kv37_to_hg19.tsv
* https://github.com/lindenb/jvarkit/blob/master/src/main/resources/chromnames/hg19_to_g1kv37.tsv
* [[VcfSampleRename]]
* [[BamRenameChromosomes]]
* http://plindenbaum.blogspot.fr/2013/07/g1kv37-vs-hg19.html


END_DOC
 */
@Program(name="vcfrenamechr",
	description="Convert the names of the chromosomes in a VCF file",
	keywords={"vcf","contig","chromosome","convert"},
	deprecatedMsg="use `bcftools annotate` with `--rename-chrs file` (does it handle BND ALT alleles ?)",
	modificationDate="20190411"
	)
public class ConvertVcfChromosomes extends Launcher {
	private static final Logger LOG = Logger.build(ConvertVcfChromosomes.class).make();
	

	@Parameter(names={"-f","--mapping","-m"},description=ContigNameConverter.OPT_DICT_OR_MAPPING_FILE_DESC,required=true)
	private Path mappingFile=null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path output= null;
	@Parameter(names={"--fatal"},description="exit with failure if a conversion fails. Default: skip the variant")
	private boolean exitOnFailure=false;

	
	
	public ConvertVcfChromosomes()
		{
		}
	
	public int doWork(final List<String> args) {
		if( this.mappingFile==null)
			{
			throw new JvarkitException.CommandLineError("undefined mapping file");
			}
		VCFIterator iterin=null;
		VariantContextWriter out = null;
		try
			{
			final ContigNameConverter customMapping = ContigNameConverter.fromPathOrOneDictionary(mappingFile);
			final Set<String> unseen=new HashSet<>();

			iterin = super.openVCFIterator(oneFileOrNull(args));
			out = this.output==null?
					VCFUtils.createVariantContextWriterToOutputStream(stdout()):
					VCFUtils.createVariantContextWriterToPath(this.output)
					;
				
			final VCFHeader header1=iterin.getHeader();
			final VCFHeader header2=new VCFHeader(
					header1.getMetaDataInInputOrder()
							.stream()
							.filter(L->!L.getKey().equals(VCFHeader.CONTIG_KEY))
							.collect(Collectors.toSet()),
						header1.getSampleNamesInOrder()
						);
			
			JVarkitVersion.getInstance().addMetaData(this, header2);
			
			if(header1.getSequenceDictionary()!=null)
				{
				header2.setSequenceDictionary(customMapping.convertDictionary(header1.getSequenceDictionary()));
				}
			out.writeHeader(header2);
			
			final Function<String, String> converter= (S)->{
				final String newName= customMapping.apply(S);
				if(StringUtils.isBlank(newName))
					{
					if(exitOnFailure) throw new IllegalArgumentException("Cannot convert contig "+S);
					if(unseen.size()<1000 && !unseen.contains(S))
						{
						LOG.warn("Cannot find contig for "+S);
						unseen.add(S);
						}
					return null;
					}
				return newName;
				};
			
			while(iterin.hasNext()) {
				final VariantContext ctx=iterin.next();
				final String contig = converter.apply(ctx.getContig());
				if(StringUtils.isBlank(contig)) continue;
					
				final VariantContextBuilder vcb= new VariantContextBuilder(ctx);

				// if it's a breaking end, loop over the ALT alleles and try to change it
				if(StructuralVariantType.BND.equals(ctx.getStructuralVariantType())) {
					boolean ok_convert = true;
					final List<Genotype> genotypes = new ArrayList<>(ctx.getGenotypes());
					final List<Allele> alleles = new ArrayList<>(ctx.getAlleles());
					for(int i=1;i< ctx.getAlleles().size();i++)
						{
						Allele alt = ctx.getAlleles().get(i);
						if(!alt.isSymbolic()) continue;
						if(!alt.isBreakpoint()) continue;
						String display = alt.getDisplayString();
						int colon = display.indexOf(':');
						if(colon==-1) continue;
						final char delim = display.contains("[")?'[':']';
						final int i1 = display.indexOf(delim);
						if(i1>colon) continue;
						final int i2 = colon+1==display.length()?-1:display.indexOf(delim,colon+1);
						if(!(i1+1 <colon && colon+1 <i2) ) {
							continue;
							}
						final String bndCtg= converter.apply(display.substring(i1+1,colon));
						if(StringUtils.isBlank(bndCtg)) {
							ok_convert=false;
							break;
							}
						Allele newAlt= Allele.create(
							display.substring(0,i1+1)+
							bndCtg + 
							display.substring(colon)
							,false);
						alleles.set(i, newAlt);
						for(int x=0;x< genotypes.size();++x) {
							final Genotype gt = genotypes.get(x);
							final List<Allele> getAlleles= new ArrayList<>(gt.getAlleles());
							getAlleles.replaceAll(A->A.equals(alt)?newAlt:A);
							genotypes.set(x, new GenotypeBuilder(gt).alleles(getAlleles).make());
							}
						}
					
					if(!ok_convert) continue;
					
					vcb.alleles(alleles);
					vcb.genotypes(genotypes);
					
				}
				
				
				vcb.chr(contig);
				out.add(vcb.make());
				}
			iterin.close();iterin=null;
			out.close();out=null;
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iterin);
			CloserUtil.close(out);
			}
		}

	public static void main(final String[] args)
		{
		new ConvertVcfChromosomes().instanceMainWithExit(args);
		}
	}
