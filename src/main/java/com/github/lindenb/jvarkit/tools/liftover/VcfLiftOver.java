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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.liftover;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
/**

## See also

picard LiftOverVcf (loads all the genome in memory...)

END_DOC
*/
@Program(name="vcfliftover",description="Lift-over a VCF file")
public class VcfLiftOver extends Launcher
	{

	private static final Logger LOG = Logger.build(VcfLiftOver.class).make();


	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	@Parameter(names={"-f","--chain"},description="LiftOver file.",required=true)
	private File liftOverFile = null;
	@Parameter(names={"-x","--failed"},description="(file.vcf) write variants failing the liftOver here. Optional.")
	private File failedFile = null;
	@Parameter(names={"-m","--minmatch"},description="lift over min-match.")
	private double userMinMatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;
	@Parameter(names={"-D","-R","-r","--reference"},description="indexed REFerence file.",required=true)
	private File faidx = null;
	@Parameter(names={"-T","--tag"},description="INFO tag")
	private String infoTag = "LIFTOVER";
	@Parameter(names={"-failtag","--failtag"},description="failed INFO tag")
	private String failedinfoTag = "LIFTOVER_FAILED";
	@Parameter(names={"-check","--check"},description="Check variant allele sequence is the same on REF")
	private boolean checkAlleleSequence = false;
	@Parameter(names={"--chainvalid"},description="Ignore LiftOver chain validation")
	private boolean ignoreLiftOverValidation=false;


	private LiftOver liftOver=null;
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	
	protected Allele revcomp(final Allele a)
		{
		if(a.isNoCall()) return a;
		if(a.isSymbolic()) return a;
		final String seq=a.getBaseString();
		final StringBuilder sb=new StringBuilder(seq.length());
		for(int i=seq.length()-1;i>=0;--i)
			{
			sb.append(AcidNucleics.complement(seq.charAt(i)));
			}
		return Allele.create(sb.toString(), a.isReference());
		}
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		VariantContextWriter failed=null;
		GenomicSequence genomicSequence = null;
		try {
			final VCFHeader incomingHeader= in.getHeader();
			if(this.failedFile!=null)
				{
				VCFHeader header2=new VCFHeader(incomingHeader);
				header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
				header2.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
				header2.addMetaDataLine(new VCFInfoHeaderLine(this.failedinfoTag,1,VCFHeaderLineType.String,"Why the liftOver failed."));

				failed= super.openVariantContextWriter(failedFile);
				failed.writeHeader(header2);
				}
			
			final VCFHeader header3=new VCFHeader(incomingHeader);
			header3.setSequenceDictionary(this.indexedFastaSequenceFile.getSequenceDictionary());
			
			
			header3.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header3.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			header3.addMetaDataLine(new VCFInfoHeaderLine(this.infoTag,1,VCFHeaderLineType.String,"Chromosome|Position before liftOver."));
			out.writeHeader(header3);
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(in.getHeader());
			while(in.hasNext())
				{
				final VariantContext ctx=progress.watch(in.next());
				
				final Interval lifted=liftOver.liftOver(
						new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd(),
						false,//negative strand
						String.join("|",ctx.getContig(),String.valueOf(ctx.getStart()),ctx.getReference().toString()))
						);
				if(lifted==null )
					{
					if(failed!=null) failed.add(new VariantContextBuilder(ctx).attribute(this.failedinfoTag, "LiftOverFailed").make());
					}
				else if(this.checkAlleleSequence && this.indexedFastaSequenceFile.getSequenceDictionary().getSequence(lifted.getContig())==null)
					{
					if(failed!=null) failed.add(new VariantContextBuilder(ctx).attribute(this.failedinfoTag, "ContigMissingInRef|"+lifted.getContig()).make());
					}
				else
					{
					boolean alleleAreValidatedVsRef=true;
					//part of the code was copied from picard/liftovervcf
					final Map<Allele, Allele> reverseComplementAlleleMap = new HashMap<>();
					final List<Allele> alleles = new ArrayList<Allele>();

	                for (final Allele oldAllele : ctx.getAlleles()) {
	                    if (lifted.isPositiveStrand() || oldAllele.isSymbolic() || oldAllele.isNoCall()) {
	                        alleles.add(oldAllele);
	                    }
	                    else {
	                        final Allele fixedAllele = Allele.create(SequenceUtil.reverseComplement(oldAllele.getBaseString()), oldAllele.isReference());
	                        alleles.add(fixedAllele);
	                        reverseComplementAlleleMap.put(oldAllele, fixedAllele);
	                        
	                        if(this.checkAlleleSequence) {
	                        	if(genomicSequence==null || !genomicSequence.getChrom().equals(lifted.getContig())) {
	                        		genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, lifted.getContig());
	                        		}
	                        	final String alleleStr = fixedAllele.getBaseString();
	                        	int x=0;
	                        	while(x<alleleStr.length() && lifted.getStart()-1+x < genomicSequence.length())
	                        		{
	                        		char refChar= genomicSequence.charAt(lifted.getStart()-1+x);
	                        		if(Character.toLowerCase(refChar)!=Character.toLowerCase(alleleStr.charAt(x)))
	                        			{
	                        			alleleAreValidatedVsRef=false;
	                        			break;
	                        			}
	                        		++x;
	                        		}
	                        	if(x!=alleleStr.length())
	                        		{
	                        		alleleAreValidatedVsRef=false;
	                        		break;
	                        		}
	                        }
	                    }
	                }
	                
	                if(!alleleAreValidatedVsRef)
	                	{
	                	if(failed!=null) failed.add(new VariantContextBuilder(ctx).attribute(this.failedinfoTag, "AlleleMismatchRef").make());
	                	continue;
	                	}

					
					final VariantContextBuilder vcb=new VariantContextBuilder(
							ctx.getSource(),
							ctx.getContig(),
							ctx.getStart(),
							ctx.getEnd(),
							alleles
							);
					vcb.id(ctx.getID());
					vcb.attributes(ctx.getAttributes());
					vcb.attribute(this.infoTag,ctx.getContig()+"|"+ctx.getStart()+"|"+ctx.getReference().getDisplayString());
					vcb.filters(ctx.getFilters());
					vcb.log10PError(ctx.getLog10PError());
					  
					final GenotypesContext genotypeContext = ctx.getGenotypes();
					final GenotypesContext fixedGenotypes = GenotypesContext.create(genotypeContext.size());
			        for ( final Genotype genotype : genotypeContext ) 
			        	{
			            final List<Allele> fixedAlleles = new ArrayList<Allele>();
			            for ( final Allele allele : genotype.getAlleles() ) {
			                final Allele fixedAllele = reverseComplementAlleleMap.containsKey(allele) ?
			                		reverseComplementAlleleMap.get(allele) : 
			                		allele;
			                fixedAlleles.add(fixedAllele);
			            	}
			            fixedGenotypes.add(new GenotypeBuilder(genotype).alleles(fixedAlleles).make());
			        	}
			        vcb.genotypes(fixedGenotypes);
				      
					
					out.add(vcb.make());
					}
				}
			if(failed!=null)
				{
				failed.close();
				failed=null;
				}
			return RETURN_OK;
		}catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(failed);	
			}
		}
	
	@Override
	public int doWork(List<String> args) {		
		if(this.liftOverFile==null)
			{
			LOG.error("LiftOver file is undefined.");
			return -1;
			}
		
		if(this.faidx==null)
			{
			LOG.error("reference file is undefined.");
			return -1;
			}
		
		try {
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(this.faidx);
			this.liftOver=new LiftOver(this.liftOverFile);
			this.liftOver.setLiftOverMinMatch(this.userMinMatch);
			if(!this.ignoreLiftOverValidation) {
				this.liftOver.validateToSequences(this.indexedFastaSequenceFile.getSequenceDictionary());
				}
			return doVcfToVcf(args,outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		
		}

	public static void main(String[] args)
		{
		new VcfLiftOver().instanceMainWithExit(args);
		}

	}
