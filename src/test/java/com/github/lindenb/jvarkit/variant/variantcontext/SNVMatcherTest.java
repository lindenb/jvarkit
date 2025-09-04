package com.github.lindenb.jvarkit.variant.variantcontext;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class SNVMatcherTest {
	@Test
	public void testID() {
		VariantContext vc1 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1)
				.id("ID1;ID2").
				alleles("A","C","G").make();
		
		VariantContext vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1)
				.id("ID3;ID2;ID4").
				alleles("A","C","G").make();
		
		
		Assert.assertTrue(SNVMatcher.id.test(vc1, vc2));
		
		vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1)
				.id("IDX;IDY").
				alleles("A","C","G").make();
		Assert.assertFalse(SNVMatcher.id.test(vc1, vc2));
		}
	
	@Test
	public void testOverlap() {
		VariantContext vc1 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","C","G").make();
		
		VariantContext vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","C","G").make();
		
		Assert.assertTrue(SNVMatcher.overlap.test(vc1, vc2));
		}
	
	@Test
	public void testChromPos() {
		VariantContext vc1 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","C","G").make();
		
		VariantContext vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("C","T","G").make();
		
		Assert.assertTrue(SNVMatcher.chrom_pos.test(vc1, vc2));
		
		vc2 = new VariantContextBuilder().
				chr("chr1").start(2).stop(2).
				alleles("C","T","G").make();
		Assert.assertFalse(SNVMatcher.chrom_pos.test(vc1, vc2));
		}
	
	@Test
	public void testChromPosRef() {
		VariantContext vc1 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","C","G").make();
		
		VariantContext vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","T","G").make();
		
		Assert.assertTrue(SNVMatcher.chrom_pos_ref.test(vc1, vc2));
		
		vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("C","T","G").make();
		Assert.assertFalse(SNVMatcher.chrom_pos_ref.test(vc1, vc2));
		}
	
	@Test
	public void testChromPosRefAnyAlt() {
		VariantContext vc1 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","C","*").make();
		
		VariantContext vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","G","C").make();
		
		Assert.assertTrue(SNVMatcher.chrom_pos_ref_any_alt.test(vc1, vc2));
		
		vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","G","T","*").make();
		Assert.assertFalse(SNVMatcher.chrom_pos_ref_any_alt.test(vc1, vc2));
		}
	
	@Test
	public void testChromPosRefAllAlt() {
		VariantContext vc1 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","C","T","*").make();
		
		VariantContext vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","C","T").make();
		
		Assert.assertTrue(SNVMatcher.chrom_pos_ref_all_alt.test(vc1, vc2));
		
		vc2 = new VariantContextBuilder().
				chr("chr1").start(1).stop(1).
				alleles("A","G","T","*").make();
		Assert.assertFalse(SNVMatcher.chrom_pos_ref_all_alt.test(vc1, vc2));
		}
	
}