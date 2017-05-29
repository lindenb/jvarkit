/**
 * 
 */
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;


import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

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
@Program(name="nozerovariationvcf",description="cat a whole VCF, or, if there is no variant, creates a fake one")
public class NoZeroVariationVCF extends Launcher
	{
	private static final Logger LOG = Logger.build(NoZeroVariationVCF.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","-r","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File faidx = null;

	
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		final VCFHeader header=in.getHeader();
		if(in.hasNext())
			{
			LOG.info("found a variant in VCF.");
			VCFUtils.copyHeaderAndVariantsTo(in, out);
			}
		else
			{
			LOG.info("no a variant in VCF. Creating a fake Variant");
			header.addMetaDataLine(new VCFFilterHeaderLine("FAKESNP", "Fake SNP created because vcf input was empty."));
			
			VCFFormatHeaderLine gtHeaderLine=header.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY);
			if(gtHeaderLine==null)
				{
				LOG.info("Adding GT to header");
				header.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
				}
			gtHeaderLine=header.getFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY);
			if(gtHeaderLine==null)
				{
				LOG.info("Adding GQ to header");
				header.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY, 1, VCFHeaderLineType.Integer, "Genotype Quality"));
				}
			
			out.writeHeader(header);
			SAMSequenceDictionary dict=this.indexedFastaSequenceFile.getSequenceDictionary();
			
			//choose random chrom, best is 'random' , but not 1...23,X,Y, etc...
			String chrom=dict.getSequence(0).getSequenceName();
			
			for(final SAMSequenceRecord ssr:dict.getSequences())
				{
				String ssn=ssr.getSequenceName();
				if(ssn.contains("_")) { chrom=ssn; break;}
				}
			
			for(final SAMSequenceRecord ssr:dict.getSequences())
				{
				String ssn=ssr.getSequenceName();
				if(ssn.toLowerCase().contains("random")) { chrom=ssn; break;}
				if(ssn.toLowerCase().contains("gl")) { chrom=ssn; break;}
				}
			
			GenomicSequence gseq=new GenomicSequence(this.indexedFastaSequenceFile,
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
			VariantContextBuilder vcb=new VariantContextBuilder();
			
			Allele a1=Allele.create((byte)ref,true);
			Allele a2=Allele.create((byte)alt,false);
			List<Allele> la1a2=new ArrayList<Allele>(2);
			List<Genotype> genotypes=new ArrayList<Genotype>(header.getSampleNamesInOrder().size());
			la1a2.add(a1);
			la1a2.add(a2);
			
			
			vcb.chr(gseq.getChrom());
			vcb.start(POS+1);
			vcb.stop(POS+1);
			vcb.filter("FAKESNP");
			vcb.alleles(la1a2);
			vcb.log10PError(-0.1);
			for(String sample:header.getSampleNamesInOrder())
				{
				final GenotypeBuilder gb=new GenotypeBuilder(sample);
				
				if(genotypes.isEmpty()) {
					gb.DP(1);
					gb.GQ(1);
					gb.alleles(la1a2);
					gb.noAD();
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
	@Override
	public int doWork(List<String> args) {
		
		if(faidx==null)
			{
			LOG.error("Indexed fasta file missing.");
			return -1;
			}
		try
			{
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
			
			return doVcfToVcf(args,this.outputFile);
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}	
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			this.indexedFastaSequenceFile =null;
			}
		}

	public static void main(String[] args)
		{
		new NoZeroVariationVCF().instanceMainWithExit(args);
		}	

}
