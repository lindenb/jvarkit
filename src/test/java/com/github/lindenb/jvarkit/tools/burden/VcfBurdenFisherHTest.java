package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class VcfBurdenFisherHTest extends TestUtils{
	@Test(dataProvider="all-vcf-files")
	public void test01(final String inputFile) 
		throws IOException
		{
		final VCFFileReader r0 = new VCFFileReader(new File(inputFile),false);
		final VCFHeader vcfheader = r0.getFileHeader();
		if(vcfheader.getNGenotypeSamples()<2) {
			r0.close();
			return;
		}
		final CloseableIterator<VariantContext> iter = r0.iterator();
		final File inputVcf = super.createTmpFile(".vcf");
		final VariantContextWriter w= VCFUtils.createVariantContextWriter(inputVcf);
		w.writeHeader(vcfheader);
		iter.stream().
			filter(V->V.getAlleles().size()==2).
			forEach(V->w.add(V));
		iter.close();
		w.close();
		r0.close();

		
		final File ped = super.createRandomPedigreeFromFile(inputVcf.getPath());
		if(ped==null) {
			 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
			return;
		}
		final File output = super.createTmpFile(".vcf");
		
		Assert.assertEquals(new VcfBurdenFisherH().instanceMain(
        		newCmd().add(
        		"-o",output,
        		"--pedigree",ped,
        		inputFile).make()
        	),0);
        assertIsVcf(output);
		}

}
