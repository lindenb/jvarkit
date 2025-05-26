package com.github.lindenb.jvarkit.tools.kg2bed;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class KnownGenesToBedTest{

@DataProvider(name = "src1")
public Object[][] createData1() {
	final TestSupport support = new TestSupport();
	return new Object[][]{
		{support.resource("rotavirus_rf.knowngenes.tsv.gz"),""},
		{support.resource("rotavirus_rf.knowngenes.tsv.gz"),"INTRON,UTR"}
		};
	}
@Test(dataProvider="src1")
public void test(final String kgfile,final String exclude) throws IOException {
	final TestSupport support = new TestSupport();
	try {
		final Path out =support.createTmpPath(".bed");
		Assert.assertEquals(new KnownGenesToBed().instanceMain(new String[] {
				"-o",out.toString(),
				"--hide",exclude,
				kgfile
				}),0);
		Assert.assertTrue(Files.exists(out));
		support.assertIsBed(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}
}
