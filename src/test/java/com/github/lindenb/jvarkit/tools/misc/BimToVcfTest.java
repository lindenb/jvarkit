package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class BimToVcfTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name="src1")
	public Object[][] getLines() {
		return new Object[][] {
			{"RF01 ID1 0 2 G A"},
			{"RF01 ID1 0 2 G GA"},
			{"RF01 ID1 0 2 GC G"},
			{"RF01 ID1 0 2 G -"},
			{"RF01 ID1 0 2 - -"},
			{"RF01 ID1 0 2 G G"}
		};
	}
	
	
	@Test(dataProvider="src1")
	public void test01(String line) throws IOException {
		try  {
			Path p = support.createTmpPath(".bim");
			Path out = support.createTmpPath(".vcf");
			Files.write(p, Collections.singletonList(line.replace(" ", "\t")));
			Assert.assertEquals(
					new BimToVcf().instanceMain(new String[] {
							"-o",out.toString(),
							"-R",support.resource("rotavirus_rf.fa"),
							p.toString()
						}),0);
			support.assertIsVcf(out);
		} finally {
			support.removeTmpFiles();
		}
		
	}
	
}
