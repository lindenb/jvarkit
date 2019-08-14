package com.github.lindenb.jvarkit.tools.gnomad;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


@AlsoTest(LauncherTest.class)
public class VcfGnomadSVTest {
	private final TestSupport support = new TestSupport();

	
@DataProvider(name="src01")
public Object[][] testData01() {
	return new Object[][] {
		{support.resource("test_vcf01.vcf")}
		
	};
}

@Test(dataProvider="src01")
public void test01(final String vcfpath) throws IOException {
	try {
		Path gnomad_sv = Paths.get("/home/lindenb/data/gnomad_v2_sv.sites.vcf.gz");
		if(!Files.exists(gnomad_sv)) return;
		final Path vcfOut = support.createTmpPath(".vcf");
		
		Assert.assertEquals(new VcfGnomad().instanceMain(new String[]{
				"-o",vcfOut.toString(),
				"-g",gnomad_sv.toString(),
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
