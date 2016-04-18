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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 *     
 * VcfBurdenRscriptV
 * @author lindenb
 *
 */
public class VcfBurdenRscriptV
	extends AbstractVcfBurdenRscriptV
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfBurdenRscriptV.class);
	
	
	public VcfBurdenRscriptV()
		{
		}
	 
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		/* 
		 mail matile April 11: requires calling R for SKAT: needs
		 matrix of samples and genotypes */
		PrintWriter pw = null;
		VcfIterator in = null;
		try {
			in = super.openVcfIterator(inputName);
			
			final VCFHeader header=in.getHeader();
			
			
			final Set<String> samplesNames= new HashSet<>(header.getSampleNamesInOrder());
			final Pedigree pedigree = Pedigree.readPedigree(header);
			if(pedigree.isEmpty()) {
				return wrapException("No pedigree found in header "+inputName+". use VcfInjectPedigree to add it");
				}
			
			if(!pedigree.verifyPersonsHaveUniqueNames()) {
				return wrapException("I can't use this pedigree in VCF because two samples have the same ID in  "+inputName);
			}
			
			//cleanup individuals
			final Set<Pedigree.Person> samples = new TreeSet<>(pedigree.getPersons());
			Iterator<Pedigree.Person> iter= samples.iterator();
			while(iter.hasNext())
			{
				final Pedigree.Person person = iter.next();
				if(!(samplesNames.contains(person.getId()) && (person.isAffected() || person.isUnaffected()))) {
					LOG.warn("Ignoring "+person+" because not in VCF header or status is unknown");
					iter.remove();
				}
			}

			pw = super.openFileOrStdoutAsPrintWriter();

			
			boolean first=true;
			pw.print("# samples ( 0: unaffected 1:affected)");
			pw.print("population <- c( ");
			for(final Pedigree.Person person : samples) {
				if(!first) pw.print(",");
				pw.print(person.isUnaffected()?0:1);
				first=false;
			}
			
			pw.println(")");
			
			List<Double> listOfMafs= new ArrayList<>();
			first=true;
			pw.println();
			pw.println("# genotypes as vector. Should be a multiple of length(samples).");
			pw.println("# 0 is homref (0/0), 1 is het (0/1), 2 is homvar (1/2)");
			pw.println("# if the variant contains another ALT allele: (0/2) and (2/2) are considered 0 (homref)");
			pw.println("genotypes <- c(");
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			while(in.hasNext()) {
				final VariantContext ctx = progess.watch(in.next());
				if(ctx.isFiltered() && !super.acceptFiltered) continue;
				final int n_alts = ctx.getAlternateAlleles().size();
				if( n_alts == 0) {
					LOG.warn("ignoring variant without ALT allele.");
					continue;
				}
				if( n_alts > 1) {
					LOG.warn("variant with more than one ALT. Using getAltAlleleWithHighestAlleleCount.");
				}
				final boolean is_chrom_X = (ctx.getContig().equals("X") ||  ctx.getContig().equals("chrX"));
				final Allele observed_alt = ctx.getAltAlleleWithHighestAlleleCount();
				double count_total = 0.0;
				double count_alt = 0.0;
				for (final Pedigree.Person person : samples) {
						final Genotype genotype = ctx.getGenotype(person.getId());
						if(genotype==null) {
							pw.close();pw=null;
							in.close();
							throw new IllegalStateException();
						}
						
						/* loop over alleles */
						for(final Allele a: genotype.getAlleles()) {
							/* chromosome X and male ? count half */
							if( is_chrom_X && person.isMale()) {
								count_total+=0.5;
								}
							else
								{
								count_total+=1.0;
								}
							if(a.equals(observed_alt))
								{
								count_alt++;
								}
							}
						
						
						
						if(!first) pw.print(",");
						if (genotype.isHomRef()) {
							pw.print('0');
						} else if (genotype.isHomVar() && genotype.getAlleles().contains(observed_alt)) {
							pw.print('2');
						} else if (genotype.isHet() && 
								genotype.getAlleles().contains(observed_alt) &&
								genotype.getAlleles().contains(ctx.getReference())) {
							pw.print('1');
						}
						/* we treat 0/2 has hom-ref */
						else if (genotype.isHet() && 
								!genotype.getAlleles().contains(observed_alt) &&
								genotype.getAlleles().contains(ctx.getReference())) {
							LOG.warn("Treating "+genotype+" as hom-ref (0) alt="+observed_alt);
							pw.print('0');
							
						} 
						/* we treat 2/2 has hom-ref */
						else if (genotype.isHomVar() && !genotype.getAlleles().contains(observed_alt)) {
							LOG.warn("Treating "+genotype+" as hom-ref (0) alt="+observed_alt);
							pw.print('0');
						}
						else {
							pw.print("-9");
						}
						first=false;
					}
				
				if(count_total>0)
					{
					final double maf = (double)count_alt/(double)count_total;
					listOfMafs.add(maf);
					}
				else
					{
					listOfMafs.add(null);
					}
				}
			pw.println(")");	
			first = true;
			
			pw.println();
			pw.println("# MAF for variants Should be a modulo==0 of length(genotypes).");
			pw.println("# 'NA' is used if no ALT is found in the VCF line.");
			pw.print("MAFs <- c(");
			for(final Double d: listOfMafs) {
				if(!first) pw.print(",");
				pw.print(d==null?"NA":String.valueOf(d));
				first=false;
			}
			pw.println(")");
			
			
			progess.finish();
			pw.flush();
			if(pw.checkError()) {
				return wrapException("I/O error");
			}
			pw.close();pw=null;
			
					
			
			
			LOG.info("done");
			return RETURN_OK;
			} catch(Exception err) {
				return wrapException(err);
			} finally {
				CloserUtil.close(pw);
				CloserUtil.close(in);
			}
		}
	
	public static void main(String[] args)
		{
		new VcfBurdenRscriptV().instanceMainWithExit(args);
		}
	}
