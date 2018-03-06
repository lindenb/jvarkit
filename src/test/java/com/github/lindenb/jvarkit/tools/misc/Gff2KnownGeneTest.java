package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

public class Gff2KnownGeneTest  extends TestUtils {
	
	@DataProvider(name="gff-data")
	public Object[][] getGffData() {
		return new Object[][] {
			{SRC_TEST_RESOURCE+"/gencode.v19.annotation.gff3"}
		};
	}
	
	@Test(dataProvider="gff-data")
	public void test01(final String inputFile) 
		throws IOException
		{
		final File out = super.createTmpFile(".kg");
		Assert.assertEquals(0,new Gff2KnownGene().instanceMain(new String[] {
			"-o",out.getPath(),
			inputFile
			}));
		final BufferedReader r = IOUtils.openFileForBufferedReading(out);
		Assert.assertTrue(r.lines().map(L->new KnownGene(L.split("[\t]"))).count()>0);
		r.close();
		}
	}
