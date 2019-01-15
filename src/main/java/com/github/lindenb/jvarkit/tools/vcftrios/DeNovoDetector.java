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
package com.github.lindenb.jvarkit.tools.vcftrios;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;

/** utility detecting denovo mutations */
public class DeNovoDetector {
	
public static interface DeNovoMutation
	{
	public Genotype getFatherGenotype();
	public Genotype getMotherGenotype();
	public Genotype getChildGenotype();
	public default boolean hasFatherGenotype() { return this.getFatherGenotype()!=null;}
	public default boolean hasMotherGenotype() { return this.getMotherGenotype()!=null;}
	public default boolean hasChildGenotype() { return this.getChildGenotype()!=null;}
	}

private boolean convertNoCallToHomRef = false;
private boolean fixPloidy = false;

/** shall we convert ./. to 0/0 ? can be used when using a VCF from gatk CombineVariants */
public void setConvertingNoCallToHomRef(boolean convertNoCallToHomRef) {
	this.convertNoCallToHomRef = convertNoCallToHomRef;
	}
public boolean isConvertingNoCallToHomRef() {
	return convertNoCallToHomRef;
	}

/** shall we convert haploid 'A' to 'A/*' (diploid and spanning deletion) ? can be used to fix calls in haploid regions */
public void setFixingPloidy(boolean fixPloidy) {
	this.fixPloidy = fixPloidy;
	}
public boolean isFixingPloidy() {
	return fixPloidy;
	}

public DeNovoMutation test(
		final VariantContext vc,
		final String fatherId,
		final String motherId,
		final String childId
		) {
	return test(vc,
			vc.getGenotype(fatherId),
			vc.getGenotype(motherId),
			vc.getGenotype(childId)
			);
	}

/**
 * 
 * @param vc
 * @param fatherId
 * @param motherId
 * @param childId
 * @return the DeNovoMutation or null
 * 
 * return null is the child genotype is null,
 * return null is the father AND mother are null,
 * return null is any genotype has ploidy != getPloidy()
 * 
 */
public DeNovoMutation test(
		final VariantContext vc,
		final Genotype fatherGt,
		final Genotype motherGt,
		final Genotype childGt
		) {
	final DeNovoMutationImpl mut = new DeNovoMutationImpl();
	mut.gChildOriginal = childGt;
	mut.gChild = convertGT(vc,mut.gChildOriginal);
	if(mut.gChild==null) return null;
	
	
	mut.gFatherOriginal = fatherGt;
	mut.gFather = convertGT(vc,mut.gFatherOriginal);
		
	mut.gMotherOriginal = motherGt;
	mut.gMother = convertGT(vc,mut.gMotherOriginal);

	return _test(mut);
	}

/** using this syntax (without 'vc'), we cannot fix genotype to HomRef, or fix the ploidy */
public DeNovoMutation test(final Genotype gFather,final  Genotype gMother,final  Genotype gChild) {
	final DeNovoMutationImpl mut = new DeNovoMutationImpl();

	mut.gChildOriginal = gChild;
	mut.gChild = gChild;
	if(mut.gChild==null) return null;
	
	mut.gFatherOriginal = gFather;
	mut.gFather = gFather;
	
	mut.gMotherOriginal = gMother;
	mut.gMother = gMother;
	return _test(mut);
	}



private DeNovoMutation _test(final DeNovoMutationImpl mut) {
	if(mut.gChild==null) return null;
	if(!testPloidy(mut.gChild)) return null;
	if(mut.gFather!=null && !testPloidy(mut.gFather)) return null;
	if(mut.gMother!=null && !testPloidy(mut.gMother)) return null;

	
	if(isNotCalled(mut.gFather) && isNotCalled(mut.gMother)) {
		return null;
		}
	else if(isNotCalled(mut.gFather) && isCalled(mut.gMother))
		{
		if(isDeNovoDuo(mut.gMother,mut.gChild)) return mut;
		}
	else if(isCalled(mut.gFather) && isNotCalled(mut.gMother))
		{
		if(isDeNovoDuo(mut.gFather,mut.gChild)) return mut;
		}
	else 
		{
		if(isDeNovoTrio(mut.gFather,mut.gMother,mut.gChild)) {
			return mut;
			}
		}
	return null;
	}

private boolean isDeNovoDuo(final Genotype parent,final Genotype child) {
	if(child.isNoCall() && parent.isCalled()) {
		return true;
	}
	
	return	isDeNovoDuo(
			child.getAllele(0),child.getAllele(1),
			parent.getAllele(0),parent.getAllele(1)
			);
	}

private boolean isDeNovoDuo(
		final Allele child1,final Allele child2,
		final Allele parent1,final Allele parent2
		)
	{
	return	 !(
			child1.equals(parent1) ||
			child1.equals(parent2) ||
			child2.equals(parent1) ||
			child2.equals(parent2)
			)
			;
	}

private boolean isDeNovoTrio(final Genotype father,final Genotype mother,final Genotype child)
	{
	if(child.isNoCall() && father.isCalled() && mother.isCalled()) return true;
	return	isDeNovoTrio(
			child,
			father.getAlleles(),
			mother.getAlleles()
			);
	}

private boolean isDeNovoTrio(
		final Genotype gChild,
		final List<Allele> fathers,
		final List<Allele> mothers
		)
	{
	for(int f=0;f< fathers.size();++f)
		{
		for(int m=0;m< mothers.size();++m)
			{
			final Genotype gt=GenotypeBuilder.create(
					gChild.getSampleName(),
					Arrays.asList(fathers.get(f),mothers.get(m))
					);
			if(gt.sameGenotype(gChild,true)) {
				

				return false;
				}
			}
		}
	return true;
	}

/** !isNotCalled(gt) */
private boolean isCalled(final Genotype gt) {
	return !isNotCalled(gt);
}
/** gt ==null || !gt.isCalled() */
private boolean isNotCalled(final Genotype gt) {
	return gt==null || !gt.isCalled();
}

private boolean testPloidy(final Genotype g) {
	return g!=null && g.getPloidy() == this.getPloidy();
	}

public int getPloidy() {
	return 2;
	}

private Genotype convertGT(final VariantContext vc,final Genotype gt) {
	return fixNoCall(vc,fixPloidy(vc,gt));
	}

private Genotype fixPloidy(final VariantContext vc,final Genotype gt) {
	if(!isFixingPloidy()) return gt;
	if(gt==null) return gt;
	if(gt.getPloidy()!=1) return gt;
	final List<Allele> alleles = Arrays.asList(gt.getAllele(0),Allele.SPAN_DEL);
	return new GenotypeBuilder(gt.getSampleName(), alleles).phased(false).make();
	}


private Genotype fixNoCall(final VariantContext vc,final Genotype gt) {
	if(!isConvertingNoCallToHomRef()) return gt;
	if(gt==null) return gt;
	if(!gt.isNoCall()) return gt;
	final Allele array[] = new Allele[gt.getPloidy()];
	Arrays.fill(array, vc.getReference());
	final List<Allele> homref_alleles = Arrays.asList(array);
	return new GenotypeBuilder(gt.getSampleName(), homref_alleles).make();
	}


private String gtToString(final Genotype g) {
	if(g==null) return "(null)";
	if(g.isNoCall()) return "[NO_CALL]";
	return "["+g.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(g.isPhased()?"|":"/"))+"]";
	}

private String gtToString(final Genotype original,final Genotype modified) {
	if(original==null) return gtToString(original);
	return original.getSampleName()+
			gtToString(original)+
			(
			original.sameGenotype(modified)?
			"":
			"->"+ gtToString(modified)
			)
			;
	}


private class DeNovoMutationImpl
	implements DeNovoMutation
	{
	Genotype gFatherOriginal;
	Genotype gFather;
	Genotype gMotherOriginal;
	Genotype gMother;
	Genotype gChildOriginal;
	Genotype gChild;

	@Override
	public Genotype getFatherGenotype() {
		return this.gFatherOriginal;
		}
	@Override
	public Genotype getMotherGenotype() {
		return this.gMotherOriginal;
		}
	@Override
	public Genotype getChildGenotype() {
		return this.gChildOriginal;
		}
	
	@Override
	public String toString() {
		return "Trio :"
				+ " Father "+ gtToString(gFatherOriginal, gFather)
				+ " Mother "+ gtToString(gMotherOriginal, gMother)
				+ " Child "+ gtToString(gChildOriginal, gChild)
				;
		}

	
	}

}
