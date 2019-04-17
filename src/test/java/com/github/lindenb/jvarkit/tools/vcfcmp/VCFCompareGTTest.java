package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

@AlsoTest(LauncherTest.class)
public class VCFCompareGTTest{
	private final TestSupport support =new TestSupport();

private VariantContext mute(final VariantContext ctx) {
	VariantContextBuilder vcb=new VariantContextBuilder(ctx);
	vcb.genotypes(ctx.getGenotypes().stream().
		map(G->support.random.nextDouble()<0.33?GenotypeBuilder.createMissing(G.getSampleName(),2):G).
		collect(Collectors.toList())
		);
	return vcb.make();
	}

private Path mute(final Path in) throws IOException {	
	Path outVcf = support.createTmpPath(".vcf");
	final VariantContextWriter w = VCFUtils.createVariantContextWriterToPath(outVcf);
	final VCFFileReader r = new VCFFileReader(in,true);
	w.writeHeader(r.getFileHeader());
	final CloseableIterator<VariantContext> iter = r.iterator();
	while(iter.hasNext())
		{
		final VariantContext ctx = iter.next();
		if(support.random.nextDouble()<0.1) continue;//ignore some variants
		w.add(mute(ctx));		
		}
	w.close();
	iter.close();
	r.close();
	support.assertIsVcf(outVcf);
	return outVcf;
	}
	
@Test(dataProvider = "all-vcf-files")
public void test01(final String vcfpath) throws Exception
	{
	try {
	Path outVcf = support.createTmpPath(".vcf");
	Path vcfIn = Paths.get(vcfpath);
	Path mute1 = mute(vcfIn);
	Path mute2 = mute(vcfIn);
	Assert.assertEquals(new  VCFCompareGT().instanceMain(new String[] {
			"-o",outVcf.toString(),
			mute1.toString(),mute2.toString()
		}),0);
	support.assertIsVcf(outVcf);
	} finally {
		support.removeTmpFiles();
		}
	}
	
	
}
