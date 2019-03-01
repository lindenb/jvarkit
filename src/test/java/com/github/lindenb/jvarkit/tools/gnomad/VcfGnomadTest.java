package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfGnomadTest {
	private final TestSupport support = new TestSupport();

	
@DataProvider(name="src01")
public Object[][] testData01() {
	return new Object[][] {
		{support.resource("test_vcf01.vcf")}
		
	};
}

private Path createManifest() throws IOException {
	final Path manifestFile = support.createTmpPath(".mft");
	final PrintWriter pw = new PrintWriter(Files.newOutputStream(manifestFile));
	pw.println("exome\t*\t"+support.resource("gnomad.exomes.r2.0.1.sites.vcf.gz"));
	pw.println("genome\t*\t"+support.resource("gnomad.genomes.r2.0.1.sites.1.vcf.gz"));
	pw.flush();
	pw.close();
	support.assertTsvTableIsConsitent(manifestFile, null);
	return manifestFile;
}

@Test(dataProvider="src01")
public void test01(final String vcfpath) throws IOException {
	try {
		final Path mFile = createManifest();
		final Path vcfOut = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfGnomad().instanceMain(new String[]{
				"-o",vcfOut.toString(),
				"-m",mFile.toString(),
				"-ac",
				"--gnomadFilter","MYF111",
				vcfpath
				}),0);
		support.assertIsVcf(vcfOut);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}

@Test(dataProvider="src01")
public void testStreaming(final String vcfpath) throws IOException {
	try {
		final Path mFile = createManifest();
		final Path vcfOut = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfGnomad().instanceMain(new String[] {
				"-o",vcfOut.toString(),
				"-m",mFile.toString(),
				"--streaming",
				vcfpath
			}),0);
		support.assertIsVcf(vcfOut);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}

@Test(dataProvider="src01")
public void testFilterGenotypes(final String vcfpath) throws IOException {
	try {
		final Path mFile = createManifest();
		final Path vcfOut = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfGnomad().instanceMain(new String[] {
			"-o",vcfOut.toString(),
			"-m",mFile.toString(),
			"-ac",
			"--gtfilter","MYF111",
			vcfpath
			}),0);
		support.assertIsVcf(vcfOut);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}

}
