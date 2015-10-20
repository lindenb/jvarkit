/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
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

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * @author lindenb
 *
 */
public class NoZeroVariationVCF extends AbstractNoZeroVariationVCF
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(NoZeroVariationVCF.class);

	
	 @Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  class MyCommand extends AbstractNoZeroVariationVCF.AbstractNoZeroVariationVCFCommand
	 	{

	
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(super.faidx==null)
			{
			return wrapException("undefined ref sequence");
			}
		try
			{
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(faidx);
			return doVcfToVcf(inputName);
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			this.indexedFastaSequenceFile=null;
			}
		}
	
	@Override
	protected Collection<Throwable> doVcfToVcf(final String inputName,VcfIterator in, VariantContextWriter out)
			throws IOException
			{
			VCFHeader header=in.getHeader();
			if(in.hasNext())
				{
				LOG.info("found a variant in VCF.");
				out.writeHeader(header);
				out.add(in.next());
				while(in.hasNext())
					{
					out.add(in.next());
					}
				}
			else
				{
				LOG.info("no a variant in VCF. Creating a fake Variant");
				header.addMetaDataLine(new VCFFilterHeaderLine("FAKESNP", "Fake SNP created because vcf input was empty. See "+getOnlineDocUrl()));
				
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
				
				for(SAMSequenceRecord ssr:dict.getSequences())
					{
					String ssn=ssr.getSequenceName();
					if(ssn.contains("_")) { chrom=ssn; break;}
					}
				
				for(SAMSequenceRecord ssr:dict.getSequences())
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
					GenotypeBuilder gb=new GenotypeBuilder(sample);
					gb.DP(1);
					gb.GQ(1);
					gb.alleles(la1a2);
					gb.noAD();
					gb.noPL();
					genotypes.add(gb.make());
					}
				vcb.genotypes(genotypes);
				vcb.noID();
				out.add(vcb.make());
				}
			return RETURN_OK;
			}
	
	 	}
	
	public static void main(String[] args)
		{
		new NoZeroVariationVCF().instanceMainWithExit(args);
		}	

}
