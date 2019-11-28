package com.github.lindenb.jvarkit.variant.variantcontext;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;

public class AttributeCleanerTest {
	private VariantContextBuilder makeCtx() {
		final VariantContextBuilder vcb = new VariantContextBuilder();
		vcb.chr("chr1");
		vcb.start(100);
		vcb.stop(100);
		vcb.id("rs25");
		vcb.filter("F1");
		vcb.filter("F2");
		vcb.filter("F3");
		vcb.log10PError(-100);
		List<Allele> alleles = Arrays.asList(Allele.REF_A,Allele.ALT_C);
		vcb.alleles(alleles);
		vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,1);
		vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,0.1);
		List<Genotype> L=new ArrayList<>();
		for(int i=0;i< 10;i++) {
			final GenotypeBuilder gb=new GenotypeBuilder("S"+i,alleles);
			gb.DP(10);
			gb.attribute("X1",0);
			L.add(gb.make());
			}
		vcb.genotypes(L);
		return vcb;
		}
	@Test
	public void testID() {
		Assert.assertFalse(AttributeCleaner.compile("ID").apply(makeCtx().make()).hasID());
	}
	@Test
	public void testQUAL() {
		Assert.assertFalse(AttributeCleaner.compile("QUAL").apply(makeCtx().make()).hasLog10PError());
		}
	
	@Test
	public void testAllinfo() {
		AttributeCleaner c=AttributeCleaner.compile("INFO");
		final VariantContext ctx=c.apply(makeCtx().make());
		Assert.assertFalse(ctx.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY));
		Assert.assertFalse(ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
		Assert.assertTrue(ctx.getAttributes().keySet().isEmpty());
		}
	@Test
	public void testAllFilters() {
		AttributeCleaner c=AttributeCleaner.compile("FILTER");
		final VariantContext ctx=c.apply(makeCtx().make());
		Assert.assertFalse(ctx.getFilters().contains("F1"));
		Assert.assertFalse(ctx.getFilters().contains("F2"));
		Assert.assertTrue(ctx.getFilters().isEmpty());
		}
	@Test
	public void testAllFormat() {
		AttributeCleaner c=AttributeCleaner.compile("FORMAT");
		final VariantContext ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.getNSamples()>0);
		Assert.assertTrue(ctx.getGenotypes().stream().allMatch(G->!G.hasDP()));
		Assert.assertTrue(ctx.getGenotypes().stream().allMatch(G->G.getExtendedAttributes().isEmpty()));
		}
	
	@Test
	public void testInfo() {
		AttributeCleaner c=AttributeCleaner.compile("INFO/AN");
		VariantContext ctx=c.apply(makeCtx().make());
		Assert.assertFalse(ctx.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY));
		Assert.assertTrue(ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
		
		c=AttributeCleaner.compile("^INFO/AN");
		ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY));
		Assert.assertFalse(ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));

		
		c=AttributeCleaner.compile("^INFO/AN,INFO/AF");
		ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY));
		Assert.assertTrue(ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
		}
	@Test
	public void testFilter() {
		AttributeCleaner c=AttributeCleaner.compile("FILTER/F1");
		VariantContext ctx=c.apply(makeCtx().make());
		Assert.assertFalse(ctx.getFilters().contains("F1"));
		Assert.assertTrue(ctx.getFilters().contains("F2"));
		
		c=AttributeCleaner.compile("^FILTER/F1");
		ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.getFilters().contains("F1"));
		Assert.assertFalse(ctx.getFilters().contains("F2"));

		
		c=AttributeCleaner.compile("FILTER/F1,FILTER/F2");
		ctx=c.apply(makeCtx().make());
		Assert.assertFalse(ctx.getFilters().contains("F1"));
		Assert.assertFalse(ctx.getFilters().contains("F2"));
		
		c=AttributeCleaner.compile("^FILTER/F1,FILTER/F2");
		ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.getFilters().contains("F1"));
		Assert.assertTrue(ctx.getFilters().contains("F2"));
		}
	
	@Test
	public void testFormat() {
		AttributeCleaner c=AttributeCleaner.compile("FORMAT/DP");
		VariantContext ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.getGenotypes().stream().allMatch(G->!G.hasDP()));
		
		c=AttributeCleaner.compile("FORMAT/X1");
		ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.getGenotypes().stream().allMatch(G->!G.hasAnyAttribute("X1")));
		
		c=AttributeCleaner.compile("FORMAT/X1,FORMAT/DP");
		ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.getGenotypes().stream().allMatch(G->!G.hasDP()));
		Assert.assertTrue(ctx.getGenotypes().stream().allMatch(G->!G.hasAnyAttribute("X1")));

		c=AttributeCleaner.compile("^FORMAT/X1,FORMAT/DP");
		ctx=c.apply(makeCtx().make());
		Assert.assertTrue(ctx.getGenotypes().stream().allMatch(G->G.hasDP()));
		Assert.assertTrue(ctx.getGenotypes().stream().allMatch(G->G.hasAnyAttribute("X1")));

		}
	
	}
