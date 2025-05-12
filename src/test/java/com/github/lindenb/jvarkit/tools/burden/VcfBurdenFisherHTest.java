package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;


public class VcfBurdenFisherHTest  {
	
	
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(
				support.
				allVcfOrBcf().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1",enabled=false)
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
		final VCFReader r0 = VCFReaderFactory.makeDefault().open(new File(inputFile),false);
		final VCFHeader vcfheader = r0.getHeader();
		if(vcfheader.getNGenotypeSamples()<2) {
			r0.close();
			return;
		}
		final CloseableIterator<VariantContext> iter = r0.iterator();
		final Path inputVcf = support.createTmpPath(".vcf");
		final VariantContextWriter w= VCFUtils.createVariantContextWriter(inputVcf.toFile());
		w.writeHeader(vcfheader);
		iter.stream().
			filter(V->V.getAlleles().size()==2).
			forEach(V->w.add(V));
		iter.close();
		w.close();
		r0.close();

		
		final Path ped = support.createRandomPedigreeFromFile(inputFile);
		if(ped==null) {
			 Reporter.getCurrentTestResult().setAttribute("warn", "No Pedigree for "+inputFile);
			return;
		}
		final Path output = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfBurdenFisherH().instanceMain(new String[] {
        	"-o",output.toString(),
        	"--pedigree",ped.toString()
			}),0);
		support.assertIsVcf(output);
		} finally {
			support.removeTmpFiles();
		}
		}

}
