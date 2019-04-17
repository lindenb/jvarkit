package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import com.github.lindenb.jvarkit.util.ucsc.KnownGeneTest;

@AlsoTest({LauncherTest.class,KnownGeneTest.class})
public class Gff2KnownGeneTest {
	private  final TestSupport support = new TestSupport();
	
	@DataProvider(name="gff-data")
	public Object[][] getGffData() {
		return new Object[][] {
			{support.resource("gencode.v19.annotation.gff3")}
		};
	}
	
	@Test(dataProvider="gff-data")
	public void test01(final String inputFile) 
		throws IOException
		{
		try {
		final Path out = support.createTmpPath(".kg");
		Assert.assertEquals(new Gff2KnownGene().instanceMain(new String[] {
			"-o",out.toString(),
			inputFile
			}),0);
		final BufferedReader r = IOUtils.openPathForBufferedReading(out);
		Assert.assertTrue(r.lines().map(L->new KnownGene(L.split("[\t]"))).count()>0);
		r.close();
		} finally {
			support.removeTmpFiles();
		}
		}
	}
