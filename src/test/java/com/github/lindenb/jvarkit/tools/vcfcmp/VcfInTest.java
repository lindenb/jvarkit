package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.misc.VCFShuffle;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

public class VcfInTest extends TestUtils {
	
final Allele ZORG_ALLELE = Allele.create("<ZORG>", false); 	
	
@Test(dataProvider="all-indexed-vcf-files")
public void testVcfIn(final String vcfpath) throws Exception
	{
	File vcf = makeVcfIn(vcfpath,"");
	Assert.assertTrue(super.variantStream(vcf).noneMatch(V->V.getAlleles().contains(ZORG_ALLELE)));
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfInAllAlt(final String vcfpath) throws Exception
	{
	File vcf = makeVcfIn(vcfpath,"--allalt");
	Assert.assertTrue(super.variantStream(vcf).noneMatch(V->V.getAlleles().contains(ZORG_ALLELE)));
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfInIndexed(final String vcfpath) throws Exception
	{
	File vcf = makeVcfIn(vcfpath,"--tabix --allalt");
	
	Assert.assertTrue(super.variantStream(vcf).noneMatch(V->V.getAlleles().contains(ZORG_ALLELE)));
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfInFilterIn(final String vcfpath) throws Exception
	{
	File vcf = makeVcfIn(vcfpath,"--filterin MYFILTER");
	Assert.assertTrue(super.variantStream(vcf).noneMatch(V->
		V.getFilters().contains("MYFILTER") && V.getAlleles().contains(ZORG_ALLELE)
		));
	Assert.assertTrue(super.variantStream(vcf).noneMatch(V->
		!V.getFilters().contains("MYFILTER") && !V.getAlleles().contains(ZORG_ALLELE)
		));
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfInFilterOut(final String vcfpath) throws Exception
	{
	File vcf = makeVcfIn(vcfpath,"--filterout MYFILTER");
	Assert.assertTrue(super.variantStream(vcf).
			filter(V->V.getFilters().contains("MYFILTER")).
			allMatch(V->V.getAlleles().contains(ZORG_ALLELE))
		);
	Assert.assertTrue(super.variantStream(vcf).
			filter(V->!V.getFilters().contains("MYFILTER")).
			noneMatch(V->V.getAlleles().contains(ZORG_ALLELE))
		);
	}


private File makeVcfIn(final String vcfpath,String other_args) throws Exception
	{
	File vcfDbIn = new File(vcfpath);
	File tmpIn = createTmpFile(".vcf");
	File outVcf = createTmpFile(".vcf");
	
	final VariantContextWriter w = VCFUtils.createVariantContextWriter(tmpIn);
	final VCFFileReader r = new VCFFileReader(vcfDbIn,true);
	
	w.writeHeader(r.getFileHeader());
	int n=0;
	final CloseableIterator<VariantContext> iter = r.iterator();
	
	while(iter.hasNext())
		{
		VariantContext ctx = iter.next();
		if(n++%2==0) continue;//ignore some variants
		w.add(ctx);
		w.add(ctx); // add it twice
		
		VariantContextBuilder vcb=new VariantContextBuilder(vcfpath, ctx.getContig(),ctx.getStart(),ctx.getEnd(), 
			Arrays.asList(ctx.getReference(),ZORG_ALLELE)	
			);
		vcb.genotypes(ctx.getGenotypes().stream().
				map(G->new GenotypeBuilder(G.getSampleName(),Arrays.asList(ZORG_ALLELE,ZORG_ALLELE)).make()
						).collect(Collectors.toList())
				);
		w.add(vcb.make());
		}
	w.close();
	iter.close();
	r.close();
	Assert.assertEquals(new VcfIn().instanceMain(newCmd().add(
			"-o",outVcf.getPath()).
			split(other_args).
			add(
			vcfDbIn,
			vcfpath
			).make()
		),0);
	assertIsVcf(outVcf);
	return outVcf;
	}
}
