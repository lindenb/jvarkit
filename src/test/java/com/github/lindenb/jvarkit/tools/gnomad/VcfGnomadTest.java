package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfGnomadTest extends TestUtils {
@DataProvider(name="src01")
public Object[][] testData01() {
	return new Object[][] {
		{SRC_TEST_RESOURCE+"/test_vcf01.vcf"}
		
	};
}

private File createManifest() throws IOException {
	final File manifestFile = super.createTmpFile(".mft");
	final PrintWriter pw = new PrintWriter(manifestFile);
	pw.println("exome\t*\t"+SRC_TEST_RESOURCE+"/gnomad.exomes.r2.0.1.sites.vcf.gz");
	pw.println("genome\t*\t"+SRC_TEST_RESOURCE+"/gnomad.genomes.r2.0.1.sites.1.vcf.gz");
	pw.flush();
	pw.close();
	assertTsvTableIsConsitent(manifestFile, null);
	return manifestFile;
}

@Test(dataProvider="src01")
public void test01(final String vcfpath) throws IOException {
	final File mFile = createManifest();
	final File vcfOut = super.createTmpFile(".vcf");
	
	Assert.assertEquals(new VcfGnomad().instanceMain(newCmd().
			add("-o",vcfOut.getPath()).
			add("-m",mFile.getPath()).
			add("-ac").
			add("--gnomadFilter","MYF111").
			add(vcfpath).make()
			),0);
	assertIsVcf(vcfOut);
	}

@Test(dataProvider="src01")
public void testStreaming(final String vcfpath) throws IOException {
	final File mFile = createManifest();
	final File vcfOut = super.createTmpFile(".vcf");
	
	Assert.assertEquals(new VcfGnomad().instanceMain(newCmd().
			add("-o",vcfOut.getPath()).
			add("-m",mFile.getPath()).
			add("--streaming").
			add(vcfpath).make()
			),0);
	assertIsVcf(vcfOut);
	}

}
