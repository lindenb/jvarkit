package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.File;
import java.io.IOException;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

public class VCFCompareGTTest  extends TestUtils {

private VariantContext mute(final VariantContext ctx) {
	VariantContextBuilder vcb=new VariantContextBuilder(ctx);
	vcb.genotypes(ctx.getGenotypes().stream().
		map(G->random.nextDouble()<0.33?GenotypeBuilder.createMissing(G.getSampleName(),2):G).
		collect(Collectors.toList())
		);
	return vcb.make();
	}

private File mute(final File in) throws IOException {	
	File outVcf = super.createTmpFile(".vcf");
	final VariantContextWriter w = VCFUtils.createVariantContextWriter(outVcf);
	final VCFFileReader r = new VCFFileReader(in,true);
	w.writeHeader(r.getFileHeader());
	final CloseableIterator<VariantContext> iter = r.iterator();
	while(iter.hasNext())
		{
		final VariantContext ctx = iter.next();
		if(random.nextDouble()<0.1) continue;//ignore some variants
		w.add(mute(ctx));		
		}
	w.close();
	iter.close();
	r.close();
	assertIsVcf(outVcf);
	return outVcf;
	}
	
@Test(dataProvider = "all-vcf-files")
public void test01(final String vcfpath) throws Exception
	{
	File outVcf = super.createTmpFile(".vcf");
	File vcfIn = new File(vcfpath);
	File mute1 = mute(vcfIn);
	File mute2 = mute(vcfIn);
	Assert.assertEquals(new  VCFCompareGT().instanceMain(newCmd().add(
			"-o",outVcf.getPath(),
			mute1,mute2
			).make()
		),0);
	assertIsVcf(outVcf);
	}
	
	
}
