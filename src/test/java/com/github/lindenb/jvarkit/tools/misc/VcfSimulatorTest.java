package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class VcfSimulatorTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][] {
			{support.resource("toy.fa")},
			{support.resource("rotavirus_rf.fa")}
			};
	}

@Test(dataProvider="src1")	
public void testRef(final String ref) throws IOException{
	try {
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(0,new VcfSimulator().instanceMain(new String[] {
			"-R",ref,
			"-o",out.toString(),
			}));
		support.assertIsVcf(out);
		} 
	finally {
		support.removeTmpFiles();
		}
	}
}
