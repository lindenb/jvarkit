package com.github.lindenb.jvarkit.tools.gff2kg;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Gff2KnownGeneTest {
	
	@DataProvider(name="gff-data")
	public Object[][] getGffData() {
		final TestSupport support = new TestSupport();
		return new Object[][] {
			{support.resource("gencode.v19.annotation.gff3")}
		};
	}
	
	@Test(dataProvider="gff-data")
	public void test01(final String inputFile) 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
		final Path out = support.createTmpPath(".kg");
		Assert.assertEquals(new Gff2KnownGene().instanceMain(new String[] {
			"-o",out.toString(),
			inputFile
			}),0);
		try(BufferedReader r = IOUtils.openPathForBufferedReading(out)) {
			//Assert.assertTrue(r.lines().map(L->new KnownGene(L.split("[\t]"))).count()>0);
			}
		} finally {
			support.removeTmpFiles();
		}
		}
	}
