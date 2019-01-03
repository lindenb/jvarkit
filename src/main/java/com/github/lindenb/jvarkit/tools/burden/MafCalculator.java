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


History:

*/
package com.github.lindenb.jvarkit.tools.burden;


import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class MafCalculator {
	public static final double NODATA=-1.0;
	private double count_total = 0;
	private double count_alt = 0;
	private boolean is_chrom_sexual;
	private final Allele observed_alt;
	private boolean no_call_is_homref=false;
		
	
	public static class Factory
		{
		@Parameter(names={"-nchr","--nocallhomref"},description="Consider no call as hom-ref")
		private boolean no_call_is_homref=false;		
		@Parameter(names={"-sexchr","--sexualchromosomes"},description="comma separated list of chromosomes that should be considered as sexual chromosomes/haploids")
		private String sexChromStr="chrX,chrY,X,Y";		
		private Set<String> sexChromSet=null;
		
		public void setNoCallIsHomRef(boolean no_call_is_homref) {
			this.no_call_is_homref = no_call_is_homref;
		}
		
		public boolean isNoCallIsHomRef() {
			return this.no_call_is_homref;
		}
		
		public MafCalculator create(final Allele observed_alt, final String contig)
			{
			if(sexChromSet==null) {
				sexChromSet=Arrays.asList(this.sexChromStr.split("[,]")).
						stream().filter(S->!S.isEmpty()).
						collect(Collectors.toSet());
				}
			final MafCalculator calc = new MafCalculator(observed_alt,sexChromSet.contains(contig));
			calc.no_call_is_homref = this.no_call_is_homref;
			return calc;
			}
		}
	
	public MafCalculator(final Allele observed_alt,final boolean is_chrom_sexual) {
		this.is_chrom_sexual = is_chrom_sexual;
		this.observed_alt = observed_alt;
		if(this.observed_alt==null) throw new IllegalArgumentException(
				"null allele"
				);
		if(!this.observed_alt.isNonReference()) throw new IllegalArgumentException(
				"Not an ALT allele: "+observed_alt
				);
		if(this.observed_alt.isNoCall()) throw new IllegalArgumentException(
				"No-call ALT allele: "+observed_alt
				);
	}
	
	public MafCalculator(final Allele observed_alt,final String contig) {
		this(observed_alt,contig.equalsIgnoreCase("chrX") || contig.equalsIgnoreCase("X") ||
			              contig.equalsIgnoreCase("chrY") || contig.equalsIgnoreCase("Y")
		);
	}
	
	/** For Matile 2016/11/29 */
	public void setNoCallIsHomRef(boolean no_call_is_homref) {
		this.no_call_is_homref = no_call_is_homref;
	}
	
	public boolean isNoCallIsHomRef() {
		return this.no_call_is_homref;
	}
	
	public double getCountTotal() {
		return this.count_total;
	}
	public double getCountAlt() {
		return this.count_alt;
	}
	
	/** count total == 0 */
	public boolean isEmpty(){
		return getCountTotal()==0;
	}
	
	public double getMaf() {
		if(getCountTotal()==0) return NODATA;
		return (double)getCountAlt()/(double)getCountTotal();
	}
	
	@Override
	public String toString() {
		return "alt="+getCountAlt()+"/total="+getCountTotal();
		}
	
	
	
	public void add(final Genotype genotype,boolean sample_is_male) {
		/* individual is not in vcf header */
		if(genotype==null) {						
			return;
		}
		
		if(!genotype.isCalled() ) {
			if(!this.isNoCallIsHomRef()) return;
			this.count_total+=( this.is_chrom_sexual && sample_is_male?1:2);
			return;
		}
		
		/* loop over alleles */
		for(final Allele a: genotype.getAlleles())
			addAllele(a,sample_is_male);
		}
	
	private void addAllele(final Allele a,boolean sample_is_male)
		{
		/* chromosome X and male ? count half */
		if( this.is_chrom_sexual && sample_is_male) {
			this.count_total+=0.5;
			}
		else
			{
			this.count_total+=1.0;
			}
		if(a.equals(observed_alt))
			{
			this.count_alt++;
			}
		}
	
	/** calculate MAF for all genotypes of this ctx with most frequent Allele . male=false*/
	public static MafCalculator get(final VariantContext ctx) {
		final MafCalculator mc = new MafCalculator(ctx.getAltAlleleWithHighestAlleleCount(),ctx.getContig());
		for(int i=0;i< ctx.getNSamples();++i) {
				mc.add(ctx.getGenotype(i), false);
			}
		return mc;
		}
	
}
