package com.github.lindenb.jvarkit.tools.burden;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;

public class MafCalculator {
	public static final double NODATA=-1.0;
	private int count_total = 0;
	private int count_alt = 0;
	private boolean is_chrom_X;
	private final Allele observed_alt;
	public MafCalculator(final Allele observed_alt,boolean is_chrom_X) {
		this.is_chrom_X = is_chrom_X;
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
	
	public MafCalculator(final Allele observed_alt,String contig) {
		this(observed_alt,contig.equalsIgnoreCase("chrX") || contig.equalsIgnoreCase("X"));
	}
	
	public int getCountTotal() {
		return this.count_total;
	}
	public int getCountAlt() {
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
		if(genotype==null || !genotype.isCalled() ) {						
			return;
		}
		/* loop over alleles */
		for(final Allele a: genotype.getAlleles())
			addAllele(a,sample_is_male);
		}
	
	private void addAllele(final Allele a,boolean sample_is_male)
		{
		/* chromosome X and male ? count half */
		if( this.is_chrom_X && sample_is_male) {
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
}
