package com.github.lindenb.jvarkit.tools.vcftrios;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class DeNovoDetectorTest  {

	private Allele a1 = Allele.create("A", true);
	private Allele a2 = Allele.create("T", false);
	private Allele a3 = Allele.create("G", false);

private Genotype gt(final String sample,Allele...alleles) {
	return new GenotypeBuilder(sample, Arrays.asList(alleles)).make();
}

private Genotype gtFather(Allele...alleles) {
	return gt("F",alleles);
}
private Genotype gtMother(Allele...alleles) {
	return gt("M",alleles);
}

private Genotype gtChild(Allele...alleles) {
	return gt("C",alleles);
}

private VariantContext variant(Genotype f,Genotype m,Genotype c)
	{
	final List<Genotype> L = new ArrayList<>(3);
	if(f!=null) L.add(f);
	if(m!=null) L.add(m);
	if(c!=null) L.add(c);
	final Set<Allele> alleles = new HashSet<>();
	alleles.add(a1);
	alleles.addAll(L.stream().flatMap(G->G.getAlleles().stream()).
			filter(A->!A.isNoCall()).
			collect(Collectors.toList()));
	return new VariantContextBuilder().chr("1").start(1).stop(1).
			alleles(alleles).
			genotypes(L).
			make();
	}

	
@Test
public void test01() {
	DeNovoDetector d= new DeNovoDetector();
	Assert.assertNull(d.test(variant(null,null,null),"F","M","C"));
	Assert.assertNull(d.test(variant(gtFather(a1,a2),gtMother(a1,a2),null),"F","M","C") );
	Assert.assertNull(d.test(variant(gtFather(a1,a1),gtMother(a1,a1),gtChild(a2)),"F","M","C") );// haploid
	Assert.assertNull(d.test(variant(gtFather(a1),gtMother(a1,a1),gtChild(a1,a2)),"F","M","C") );// haploid
	Assert.assertNull(d.test(variant(gtFather(a1,a1),gtMother(a1),gtChild(a1,a2)),"F","M","C") );// haploid
	Assert.assertNull(d.test(variant(gtFather(a1,a2),gtMother(a1,a2),gtChild(a1,a2)),"F","M","C") );
	
	Assert.assertNull(d.test(variant(gtFather(a1,a2),gtMother(a1,a2),gtChild(a3,a3)),"F","M","NOT_CHILDREN") );
	
	Assert.assertNull(d.test(variant(gtFather(a1,a1),gtMother(a1,a1),gtChild(a1,a1)),"F","M","C") );
	Assert.assertNull(d.test(variant(gtFather(a1,a1),gtMother(a2,a2),gtChild(a1,a2)),"F","M","C") );
	Assert.assertNull(d.test(variant(gtFather(a1,a2),gtMother(a1,a2),gtChild(a2,a2)),"F","M","C") );
	Assert.assertNull(d.test(variant(gtFather(a1,a2),gtMother(a1,a2),gtChild(a1,a1)),"F","M","C") );

	}

@Test
public void test02() {
	DeNovoDetector d= new DeNovoDetector();
	

	
	Assert.assertNotNull(d.test(variant(
			null,
			gtMother(a1,a1),
			gtChild(a2,a3) // a3 is denovo , duo
			),"F","M","C") );
	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a1),
			null,
			gtChild(a2,a3) // a3 is denovo , duo
			),"F","M","C") );

	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(a1,a3) // a3 is denovo
			),"F","M","C") );
	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(a3,a3) // a3 is denovo
			),"F","M","C") );

	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(Allele.NO_CALL,Allele.NO_CALL) // could be a deletion
			),"F","M","C") );
	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a1),
			gtMother(a1,a1),
			gtChild(a1,Allele.SPAN_DEL) // could be a haploid deletion
			),"F","M","C") );
	}
@Test
public void test03() {
	DeNovoDetector d= new DeNovoDetector();
	d.setConvertingNoCallToHomRef(true);
	
	
	Assert.assertNull(d.test(variant(
			gtFather(a1,a1),
			gtMother(a1,a1),
			gtChild(Allele.NO_CALL,Allele.NO_CALL) // -> converted to a1/a1
			),"F","M","C") );
	
	Assert.assertNull(d.test(variant(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(Allele.NO_CALL,Allele.NO_CALL) // -> converted to a1/a1
			),"F","M","C") );
	
	Assert.assertNull(d.test(variant(
			gtFather(a1,a1),
			gtMother(a1,a1),
			gtChild(Allele.NO_CALL,Allele.NO_CALL) // -> converted to a1/a1
			),"F","M","C") );

	Assert.assertNotNull(d.test(variant(
			gtFather(a2,a2),
			gtMother(a2,a2),
			gtChild(Allele.NO_CALL,Allele.NO_CALL) // -> converted to a1/a1
			),"F","M","C") );

	
	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a1),
			gtMother(a1,a1),
			gtChild(Allele.SPAN_DEL,Allele.SPAN_DEL) 
			),"F","M","C") );
	}

@Test
public void test04() {
	DeNovoDetector d= new DeNovoDetector();
	d.setFixingPloidy(true);
	
	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(a1) // -> converted to a1/<SPAN>
			),"F","M","C") );
	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(a1) // -> converted to a1/<SPAN>
			),"F","M","C") );
	
	Assert.assertNotNull(d.test(variant(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(a3) // -> converted to a3/<SPAN>
			),"F","M","C") );
	}
@Test
public void test05() {
	DeNovoDetector d= new DeNovoDetector();
	
	Assert.assertNull(d.test(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(a1,a2)
			));
	
	Assert.assertNotNull(d.test(
			gtFather(a1,a2),
			gtMother(a1,a2),
			gtChild(a1,a3)
			));
	}

}
