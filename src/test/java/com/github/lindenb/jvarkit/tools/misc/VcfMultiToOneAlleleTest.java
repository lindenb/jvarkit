package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class VcfMultiToOneAlleleTest {
	private final Allele REF = Allele.create("A", true);
	private final Allele A1 = Allele.create("G", false);
	private final Allele A2 = Allele.create("C", false);
	
	private final TestSupport support = new TestSupport();

	private List<VariantContext> createVcf(final String params,final List<Genotype> genotypes) throws IOException
		{
		try {
		Set<VCFHeaderLine> metaData = new HashSet<>();
		VCFStandardHeaderLines.addStandardFormatLines(metaData, true, "GT","DP");
		VCFStandardHeaderLines.addStandardInfoLines(metaData, true, "AF","AN","AC","DP");
		final VCFHeader vcfheader = new VCFHeader(metaData,genotypes.stream().
				map(G->G.getSampleName()).
				sorted().
				collect(Collectors.toSet())
				);
		
		final VariantAttributesRecalculator calc = new VariantAttributesRecalculator();
		calc.setHeader(vcfheader);
		final Path vcfOut= support.createTmpPath(".vcf");
		final Path vcfOut2= support.createTmpPath(".vcf");
		VariantContextBuilder vcb =new VariantContextBuilder();
		vcb.chr("1");
		vcb.start(1);
		vcb.stop(1);
		vcb.alleles(genotypes.stream().flatMap(G->G.getAlleles().stream()).collect(Collectors.toSet()));
		vcb.genotypes(genotypes);
		
		final VariantContextWriter w= VCFUtils.createVariantContextWriter(vcfOut.toFile());
		w.writeHeader(vcfheader);
		w.add(calc.apply(vcb.make()));
		w.close();
		
		List<String> args = new ArrayList<>();
		args.add("-o");args.add(vcfOut2.toString());
		Arrays.asList(params.split("[ ]")).stream().filter(S->!S.isEmpty()).forEach(S->args.add(S));
		args.add(vcfOut.toString());
		
		Assert.assertEquals(new VcfMultiToOneAllele().instanceMain(args),0);
		support.assertIsVcf(vcfOut2);
		return support.variantStream(vcfOut2).
				collect(Collectors.toList());
		} finally {
			support.removeTmpFiles();
		}
		}
	@Test
	public void testNoSampleNoMultiple() throws IOException
		{
		List<Genotype> genotypes = new ArrayList<>();
		genotypes.add(new GenotypeBuilder("S1",Arrays.asList(REF,A1)).make());
		genotypes.add(new GenotypeBuilder("S2",Arrays.asList(A1,A1)).make());
		genotypes.add(new GenotypeBuilder("S3",Arrays.asList(REF,REF)).make());
		List<VariantContext> variants = createVcf("",genotypes);
		Assert.assertEquals(variants.size(), 1);
		VariantContext ctx  = variants.get(0);
		Assert.assertEquals(ctx.getGenotypes().size(),0);
		}
	
	@Test
	public void testNoMultiple() throws IOException
		{
		List<Genotype> genotypes = new ArrayList<>();
		genotypes.add(new GenotypeBuilder("S1",Arrays.asList(REF,A1)).make());
		genotypes.add(new GenotypeBuilder("S2",Arrays.asList(A1,A1)).make());
		genotypes.add(new GenotypeBuilder("S3",Arrays.asList(REF,REF)).make());
		List<VariantContext> variants = createVcf("--samples",genotypes);
		Assert.assertEquals(variants.size(), 1);
		VariantContext ctx  = variants.get(0);
		Assert.assertEquals(ctx.getGenotypes().size(),genotypes.size());
		Assert.assertTrue(ctx.getGenotype("S1").sameGenotype(genotypes.get(0)));
		Assert.assertTrue(ctx.getGenotype("S2").sameGenotype(genotypes.get(1)));
		Assert.assertTrue(ctx.getGenotype("S3").sameGenotype(genotypes.get(2)));
		}
	@Test
	public void testSimple() throws IOException
		{
		List<Genotype> genotypes = new ArrayList<>();
		genotypes.add(new GenotypeBuilder("S1",Arrays.asList(REF,A1)).make());
		genotypes.add(new GenotypeBuilder("S2",Arrays.asList(A1,A1)).make());
		genotypes.add(new GenotypeBuilder("S3",Arrays.asList(A1,A2)).make());
		List<VariantContext> variants = createVcf("",genotypes);
		Assert.assertEquals(variants.size(), 2);
		for(VariantContext v :variants) {
			Assert.assertEquals(v.getNAlleles(),2);
			Assert.assertTrue(v.getGenotypes().stream().noneMatch(G->G.isHomVar()));
			}
		}
}
