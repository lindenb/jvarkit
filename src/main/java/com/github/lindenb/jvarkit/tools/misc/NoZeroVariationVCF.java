/**
 * 
 */
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
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
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * @author lindenb
 *
 */
public class NoZeroVariationVCF extends AbstractVCFFilter2
	{
	@Override
	public String getProgramDescription() {
		return "cat a whole VCF, or, if there is no variant, creates a fake one ";
		}
	
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/NoZeroVariationVCF";
	}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -r (file)  fasta sequence file. "+getMessageBundle("reference.faidx")+" . Required.");
		super.printOptions(out);
		}
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header=in.getHeader();
		if(in.hasNext())
			{
			info("found a variant in VCF.");
			out.writeHeader(header);
			out.add(in.next());
			while(in.hasNext())
				{
				out.add(in.next());
				}
			}
		else
			{
			info("no a variant in VCF. Creating a fake Variant");
			header.addMetaDataLine(new VCFFilterHeaderLine("FAKESNP", "Fake SNP created because vcf input was empty. See "+getOnlineDocUrl()));
			
			VCFFormatHeaderLine gtHeaderLine=header.getFormatHeaderLine(VCFConstants.GENOTYPE_KEY);
			if(gtHeaderLine==null)
				{
				info("Adding GT to header");
				header.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));
				}
			gtHeaderLine=header.getFormatHeaderLine(VCFConstants.GENOTYPE_QUALITY_KEY);
			if(gtHeaderLine==null)
				{
				info("Adding GQ to header");
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
		}
	
	@Override
	public int doWork(String[] args)
		{
		String faidx=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args, "hvL:r:"))!=-1)
			{
			switch(c)
				{
				case 'h': printUsage();return 0;
				case 'v': System.out.println(getVersion());return 0;
				case 'L': getLogger().setLevel(java.util.logging.Level.parse(getopt.getOptArg()));break;
				case 'r': faidx=getopt.getOptArg();break;
				case ':': System.err.println("Missing argument for option -"+getopt.getOptOpt());return -1;
				default: System.err.println("Unknown option -"+getopt.getOptOpt());return -1;
				}
			}
		if(faidx==null)
			{
			error("Indexed fasta file missing.");
			return -1;
			}
		try
			{
			info("opening "+faidx);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(new File(faidx));
			
			int ret= doWork(getopt.getOptInd(), args);
			
			this.indexedFastaSequenceFile.close();
			return ret;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}		
		}

	public static void main(String[] args)
		{
		new NoZeroVariationVCF().instanceMainWithExit(args);
		}	

}
