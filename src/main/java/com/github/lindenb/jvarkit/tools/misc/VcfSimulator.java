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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;


/**
BEGIN_DOC

TODO


END_DOC
*/
@Program(name="vcfSimulator",
	description="Generate a VCF",
	keywords= {"vcf"})
public class VcfSimulator extends Launcher
	{
	private static Logger LOG=Logger.build(SkipXmlElements.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names="-S",description="random number seed. -1 == current time.")
	private long randomSeed = -1L;
	@Parameter(names="-R",description=INDEXED_FASTA_REFERENCE_DESCRIPTION,converter=Launcher.IndexedFastaSequenceFileConverter.class,required=true)
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	@Parameter(names= {"-n","--samples"},description="Number of samples. Random if lower than zero")
	private int numSamples=-1;
	@Parameter(names= {"-nv","--num-variants"},description="Number of variants to output. <0 == no limit.")
	private long numberOfVariants = 1000;
	
	private Set<String> samples=new HashSet<String>();
	
	
	private Random random;
	public VcfSimulator()
		{
		}
	


	@Override
	public int doWork(final List<String> args)
		{
		if(this.indexedFastaSequenceFile==null)
			{
			LOG.error("Reference is undefined");
			return -1;
			}
		if(!args.isEmpty()) {
			LOG.error("too many arguments");
			return -1;
		}
		
		if(this.randomSeed==-1L) {
			this.random = new Random(System.currentTimeMillis());
			}
		else
			{
			this.random = new Random(this.randomSeed);
			}
		if(this.numSamples<0)
			{
			this.numSamples=1+this.random.nextInt(10);
			}
		
		while(this.samples.size()<numSamples) this.samples.add("SAMPLE"+(1+this.samples.size()));
		VariantContextWriter writer=null;
		PrintStream pw = null;
		try
			{
			final Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true, "GT","DP");
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true, "AF","AN","AC","DP");
			VariantAttributesRecalculator calc= new VariantAttributesRecalculator();
			
			final VCFHeader header=new VCFHeader(
					metaData,
					this.samples);
			header.setSequenceDictionary(this.indexedFastaSequenceFile.getSequenceDictionary());
			calc.setHeader(header);
				
			pw = super.openFileOrStdoutAsPrintStream(this.outputFile);
			writer = VCFUtils.createVariantContextWriterToOutputStream(pw);
			
			writer.writeHeader(header);
			long countVariantsSoFar=0;
			for(;;)
				{
				if(pw.checkError()) break;
				if(this.numberOfVariants>=0 && countVariantsSoFar>=this.numberOfVariants) break;
				for(final SAMSequenceRecord ssr: this.indexedFastaSequenceFile.getSequenceDictionary().getSequences())
					{
					if(pw.checkError()) break;
					final GenomicSequence genomic=new GenomicSequence(this.indexedFastaSequenceFile, ssr.getSequenceName());
					for(int pos=1;pos<=ssr.getSequenceLength();++pos)
						{
						if(pw.checkError()) break;
						char REF=Character.toUpperCase(genomic.charAt(pos-1));
						if(REF=='N') continue;
						char ALT='N';
						switch(REF)
							{
							case 'A': ALT="TGC".charAt(random.nextInt(3));break;
							case 'T': ALT="AGC".charAt(random.nextInt(3));break;
							case 'G': ALT="TAC".charAt(random.nextInt(3));break;
							case 'C': ALT="TGA".charAt(random.nextInt(3));break;
							default: ALT='N';
							}
						if(ALT=='N') continue;
						final Allele refAllele=Allele.create((byte)REF,true);
						Allele altAllele=Allele.create((byte)ALT,false);
						final VariantContextBuilder cb=new VariantContextBuilder();
						cb.chr(genomic.getChrom());
						cb.start(pos);
						cb.stop(pos);
						List<Genotype> genotypes=new ArrayList<Genotype>(samples.size());
						for(String sample: samples)
							{
							final Allele a1=(random.nextBoolean()?refAllele:altAllele);
							final Allele a2=(random.nextBoolean()?refAllele:altAllele);
							GenotypeBuilder gb=new GenotypeBuilder(
									sample,
									Arrays.asList(a1,a2)
									);
							if(random.nextBoolean())
								{
								gb=new GenotypeBuilder(
										sample,
										Arrays.asList(a1,a2)
										);
								gb.DP(1+random.nextInt(50));
								
								}
							genotypes.add(gb.make());
							}
						cb.genotypes(genotypes);
						cb.alleles(Arrays.asList(refAllele,altAllele));
						writer.add(calc.apply(cb.make()));
						countVariantsSoFar++;
						}
					
					}
				}
			writer.close();
			writer = null;
			pw.flush();
			pw.close();
			pw=null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			CloserUtil.close(pw);
			CloserUtil.close(writer);
			}
		}

	public static void main(String[] args)
		{
		new VcfSimulator().instanceMainWithExit(args);
		}
	}
