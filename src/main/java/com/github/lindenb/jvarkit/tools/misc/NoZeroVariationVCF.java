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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceFileSupplier;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example

```bash
# we use grep to create an empty VCF
$ gunzip -c file.vcf.gz | \
 grep  "#" |\
 java -jar dist/nozerovariationvcf.jar -r human_g1k_v37.fasta

##fileformat=VCFv4.1
##FILTER=<ID=FAKESNP,Description="Fake SNP created because vcf input was empty. See https://github.com/lindenb/jvarkit">
(...)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
GL000207.1	1	.	C	A	1	FAKESNP	.	GT:DP:GQ	0/1:1:1
```

END_DOC
*/
@Program(name="nozerovariationvcf",
description="cat a whole VCF, or, if there is no variant, creates a fake one",
keywords="vcf"
)
public class NoZeroVariationVCF extends Launcher
	{
	private static final Logger LOG = Logger.build(NoZeroVariationVCF.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","-r","--reference"},description=ReferenceFileSupplier.OPT_DESCRIPTION)
	private ReferenceFileSupplier referenceFileSupplier = ReferenceFileSupplier.getDefaultReferenceFileSupplier();
	@Parameter(names={"-f","--filter"},description="FILTER name")
	private String filter = "FAKESNP";
	
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator in,
			final VariantContextWriter out
			) {
			IndexedFastaSequenceFile indexedFastaSequenceFile=null;
			try {
				final File faidx = this.referenceFileSupplier.getRequired();
				 indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
		
				final VCFHeader header=in.getHeader();
				if(in.hasNext())
					{
					LOG.info("found a variant in VCF.");
					VCFUtils.copyHeaderAndVariantsTo(in, out);
					}
				else
					{
					LOG.info("no a variant in VCF. Creating a fake Variant");
					final Set<VCFHeaderLine> meta = new HashSet<>();
					
					meta.add(new VCFFilterHeaderLine(this.filter, "Fake SNP created because vcf input was empty."));
					meta.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
					meta.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY,true));
					meta.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY,true));
					
					
					
					for(final VCFHeaderLine L:meta) header.addMetaDataLine(L);
					out.writeHeader(header);
					final SAMSequenceDictionary dict= indexedFastaSequenceFile.getSequenceDictionary();
					if(dict==null || dict.isEmpty())
						{
						throw new JvarkitException.FastaDictionaryMissing(faidx);
						}
					//choose random chrom, best is 'random' , but not 1...23,X,Y, etc...
					String chrom=dict.getSequence(0).getSequenceName();
					
					for(final SAMSequenceRecord ssr:dict.getSequences())
						{
						final String ssn=ssr.getSequenceName();
						if(ssn.contains("_")) { chrom=ssn; break;}
						}
					
					for(final SAMSequenceRecord ssr:dict.getSequences())
						{
						final String ssn=ssr.getSequenceName();
						if(ssn.toLowerCase().contains("random")) { chrom=ssn; break;}
						if(ssn.toLowerCase().contains("gl")) { chrom=ssn; break;}
						}
					
					GenomicSequence gseq=new GenomicSequence(indexedFastaSequenceFile,
							chrom
							);
					char ref='N';
					char alt='N';
					int POS=0;
					for(POS=0;POS< gseq.length();++POS)
						{
						ref=Character.toUpperCase(gseq.charAt(POS));
						if(ref=='N') continue;
						switch(ref)
							{
							case 'A': alt='T'; break;
							case 'T': alt='G'; break;
							case 'G': alt='C'; break;
							case 'C': alt='A'; break;
							default:break;
							}
						if(alt=='N') continue;
						break;
						}
					if(alt=='N') throw new RuntimeException("found only N");
					final VariantContextBuilder vcb=new VariantContextBuilder();
					
					final Allele a1=Allele.create((byte)ref,true);
					final Allele a2=Allele.create((byte)alt,false);
					final List<Allele> la1a2= Arrays.asList(a1,a2);
					
					final List<Genotype> genotypes=new ArrayList<Genotype>(header.getSampleNamesInOrder().size());
		
					vcb.chr(gseq.getChrom());
					vcb.start(POS+1);
					vcb.stop(POS+1);
					vcb.filter(this.filter);
					vcb.alleles(la1a2);
					vcb.log10PError(-0.1);
					for(final String sample:header.getSampleNamesInOrder())
						{
						final GenotypeBuilder gb=new GenotypeBuilder(sample);
						
						if(genotypes.isEmpty()) {
							gb.DP(0);
							gb.GQ(0);
							gb.alleles(la1a2);
							gb.AD(new int[]{0,0});
							gb.noPL();
							genotypes.add(gb.make());
							}
						else
							{
							genotypes.add(GenotypeBuilder.createMissing(sample, 2));
							}
						}
					vcb.genotypes(genotypes);
					vcb.noID();
					out.add(vcb.make());
					}
				return 0;
				}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally {
			CloserUtil.close(indexedFastaSequenceFile);
			indexedFastaSequenceFile =null;
			}
		
		}
	@Override
	public int doWork(final List<String> args) {
		return doVcfToVcf(args,this.outputFile);
		}

	public static void main(final String[] args)
		{
		new NoZeroVariationVCF().instanceMainWithExit(args);
		}	

}
