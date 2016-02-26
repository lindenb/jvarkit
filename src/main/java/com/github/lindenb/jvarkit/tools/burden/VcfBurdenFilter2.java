/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.util.Collection;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfBurdenFilter2
	extends AbstractVcfBurdenFilter2
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfBurdenFilter2.class);
	
	public VcfBurdenFilter2()
		{
		}
	 
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.casesFile==null || !super.casesFile.exists()) {
		return wrapException("Undefined Case file option -"+OPTION_CASESFILE);
		}
		if(super.controlsFile==null || !super.controlsFile.exists()) {
		return wrapException("Undefined Control file option -"+OPTION_CONTROLSFILE);
		}
		return super.initializeKnime();
	 	}

	
	
	/* public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator in,
			final VariantContextWriter out
			) throws IOException {
		final VCFHeader header=in.getHeader();
		final Set<Pedigree.Person> caseSamples = Pedigree.readPedigree(super.casesFile).getPersons();
		final Set<Pedigree.Person> controlSamples = Pedigree.readPedigree(super.controlsFile).getPersons();
		for(int pop=0;pop<2;++pop) {
			final Set<Pedigree.Person> persons = (pop==0?caseSamples:controlSamples);
			for(Pedigree.Person sample:persons) {
				if(!header.getSampleNamesInOrder().contains(sample.getId())) {
					LOG.warn("VCF header doesn't contain sample "+sample);
				}
			}
		}
		
		try {
			final VCFHeader h2=addMetaData(new VCFHeader(header));

			final VCFInfoHeaderLine mafCasInfoHeader = new VCFInfoHeaderLine(
					"BurdenF2MAFCas",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"Burden Filter F2. MAF Cases"
					);
			final VCFInfoHeaderLine mafControlsInfoHeader = new VCFInfoHeaderLine(
					"BurdenF2MAFControls",VCFHeaderLineCount.A,VCFHeaderLineType.Float,"Burden Filter F2. MAF Cases"
					);
			final VCFFilterHeaderLine filterCasHeader = new VCFFilterHeaderLine(
					mafCasInfoHeader.getID(),"MAF of cases is greater than "+super.maxMAF
					);
			final VCFFilterHeaderLine filterControlsHeader = new VCFFilterHeaderLine(
					mafControlsInfoHeader.getID(),"MAF of controls is greater than "+super.maxMAF
					);
			final VCFFilterHeaderLine filterCaseOrControlsHeader = new VCFFilterHeaderLine(
					"BurdenF2MAFCaseOrControls","MAF of (cases OR controls) is greater than "+super.maxMAF
					);

			h2.addMetaDataLine(mafCasInfoHeader);
			h2.addMetaDataLine(mafControlsInfoHeader);
			h2.addMetaDataLine(filterCasHeader);
			h2.addMetaDataLine(filterControlsHeader);
			h2.addMetaDataLine(filterCaseOrControlsHeader);
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				final VariantContext ctx = progess.watch(in.next());
				
				
				if(	ctx.getAlternateAlleles().size()!=1) {
					return wrapException("Expected only one allele per variant. Please use ManyAlleletoOne.");
					}
				final Allele observed_alt = ctx.getAlternateAllele(0);
				final boolean is_chrom_X = (ctx.getContig().equals("X") ||  ctx.getContig().equals("chrX"));

				int count_filter_set=0;
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				/* loop over two populations : 0 = case, 1=controls */
				for(int pop=0;pop<2;++pop) {
					int count_alt = 0;
					double count_total=0.0;/* total number of alleles */
					/* loop over persons in this pop */
					for(final Pedigree.Person p:(pop==0?caseSamples:controlSamples)) 	{
					/* get genotype for this individual */
					final Genotype genotype = ctx.getGenotype(p.getId());
					/* individual is not in vcf header */
					if(genotype==null) continue;
					/* genotype is not called */
					if(!genotype.isCalled()) continue;
					/* loop over alleles */
					for(final Allele a: genotype.getAlleles()) {
						/* chromosome X and male ? count half */
						if( is_chrom_X && p.getSex()==Pedigree.Sex.male) {
							count_total+=0.5;
							}
						else
							{
							count_total+=1.0;
							}
						if(a.equals(observed_alt)) count_alt++;//TODO even on X ??
						}
					}
					/* at least one genotype foud */
					if(count_total!=0)
						{
						/* get MAF */
						double maf = (double)count_alt/(double)count_total;
						/* set INFO attribute */
						vcb.attribute(
								(pop==0?mafCasInfoHeader.getID():mafControlsInfoHeader.getID()),	
								maf);
						/* set FILTER if needed */
						if(maf>super.maxMAF) {
							vcb.filter(pop==0?filterCasHeader.getID():filterControlsHeader.getID());
							count_filter_set++;
							}
						} 
					}/* end of loop over pop */
				if(count_filter_set==2) {
					vcb.filter(filterCaseOrControlsHeader.getID());
				}
				out.add(vcb.make());
				}
			progess.finish();
			return RETURN_OK;
			} catch(Exception err) {
				return wrapException(err);
			} finally {
				CloserUtil.close(in);
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenFilter2().instanceMainWithExit(args);
		}
	}
