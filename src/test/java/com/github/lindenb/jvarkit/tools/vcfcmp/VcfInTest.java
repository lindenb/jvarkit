package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFReader;


public class VcfInTest  {
private final TestSupport support =new TestSupport();
	
final Allele ZORG_ALLELE = Allele.create("<ZORG>", false); 	

@DataProvider(name = "all-indexed-vcf-files")
public Object[][] createData1() {
	return new Object[][] {
			{support.resource("test_vcf01.vcf")},
			{support.resource("manta.B00GWGD.vcf.gz")},
			{support.resource("rotavirus_rf.vcf.gz")},
			{support.resource("S5.vcf.gz")},
			{support.resource("S4.vcf.gz")},
			{support.resource("toy.vcf.gz")},
			{support.resource("rotavirus_rf.ann.vcf.gz")},
			{support.resource("manta.D000Q1R.vcf.gz")},
			{support.resource("S2.vcf.gz")},
			{support.resource("rotavirus_rf.freebayes.vcf.gz")},
			{support.resource("ExAC.r1.sites.vep.vcf.gz")},
			{support.resource("gnomad.genomes.r2.0.1.sites.1.vcf.gz")},
			{support.resource("S3.vcf.gz")},
			{support.resource("S1.vcf.gz")},
			{support.resource("manta.B00GWIU.vcf.gz")},
			{support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz")},
			{support.resource("gnomad_v2_sv.sites.vcf.gz")},
			{support.resource("manta.B00I9CJ.vcf.gz")},
			{support.resource("rotavirus_rf.unifiedgenotyper.vcf.gz")},
			{support.resource("rotavirus_rf.unifiedgenotyper.vcf.gz")}
			}; 
			
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfIn(final String vcfpath) throws Exception
	{
	Path vcf = makeVcfIn(vcfpath,"");
	Assert.assertTrue(support.variantStream(vcf).noneMatch(V->V.getAlleles().contains(ZORG_ALLELE)));
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfInAllAlt(final String vcfpath) throws Exception
	{
	Path vcf = makeVcfIn(vcfpath,"--allalt");
	Assert.assertTrue(support.variantStream(vcf).noneMatch(V->V.getAlleles().contains(ZORG_ALLELE)));
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfInIndexed(final String vcfpath) throws Exception
	{
	Path vcf = makeVcfIn(vcfpath,"--tabix --allalt");
	
	Assert.assertTrue(support.variantStream(vcf).noneMatch(V->V.getAlleles().contains(ZORG_ALLELE)));
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfInFilterIn(final String vcfpath) throws Exception
	{
	Path vcf = makeVcfIn(vcfpath,"--filterin MYFILTER");
	Assert.assertTrue(support.variantStream(vcf).noneMatch(V->
		V.getFilters().contains("MYFILTER") && V.getAlleles().contains(ZORG_ALLELE)
		));
	Assert.assertTrue(support.variantStream(vcf).noneMatch(V->
		!V.getFilters().contains("MYFILTER") && !V.getAlleles().contains(ZORG_ALLELE)
		));
	}

@Test(dataProvider="all-indexed-vcf-files")
public void testVcfInFilterOut(final String vcfpath) throws Exception
	{
	Path vcf = makeVcfIn(vcfpath,"--filterout MYFILTER");
	Assert.assertTrue(support.variantStream(vcf).
			filter(V->V.getFilters().contains("MYFILTER")).
			allMatch(V->V.getAlleles().contains(ZORG_ALLELE))
		);
	Assert.assertTrue(support.variantStream(vcf).
			filter(V->!V.getFilters().contains("MYFILTER")).
			noneMatch(V->V.getAlleles().contains(ZORG_ALLELE))
		);
	}


private Path makeVcfIn(final String vcfpath,String other_args) throws Exception
	{
	Path vcfDbIn = Paths.get(vcfpath);
	Path tmpIn = support.createTmpPath(".vcf");
	Path outVcf = support.createTmpPath(".vcf");
	
	final VariantContextWriter w = VCFUtils.createVariantContextWriterToPath(tmpIn);
	final VCFReader r = VCFReaderFactory.makeDefault().open(vcfDbIn,true);
	
	w.writeHeader(r.getHeader());
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
	
	final List<String> args= new ArrayList<>();
	args.add("-o");
	args.add(outVcf.toString());
	args.add("--database");
	args.add(vcfDbIn.toString());
	
	Arrays.stream(other_args.split("[ ]")).filter(S->!S.isEmpty()).forEach(S->args.add(S));
	args.add(vcfpath.toString());
	
	Assert.assertEquals(new VcfIn().instanceMain(args),0);
	support.assertIsVcf(outVcf);
	return outVcf;
	}

@Test
public void testMultiple() throws Exception
	{
	final Path vcfList = support.createTmpPath(".list");
	PrintWriter pw = new PrintWriter(Files.newBufferedWriter(vcfList));
	for(int i=1;i<6;i++)
		{
		pw.println(support.resource("S"+i+".vcf.gz"));
		}
	pw.flush();
	pw.close();
	
	
	final Path outVcf = support.createTmpPath(".vcf");
	Assert.assertEquals(new VcfIn().instanceMain(new String[] {
			"-o",outVcf.toString(),
			"--database",vcfList.toString(),
			"--tabix",
			support.resource("S1.vcf.gz")
			}),0);
	support.assertIsVcf(outVcf);
	}


}
