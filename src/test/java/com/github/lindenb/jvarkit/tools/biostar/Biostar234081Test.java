package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Biostar234081Test  {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return (Object[][])support.
				allSamOrBams().
				map(F->new Object[] {F}).
				toArray()
				;
		}
	
	@Test(dataProvider="src1")
	public void test1(final String samFile) throws IOException {
		try {
		final Path out = this.support.createTmpPath(".bam");
		Assert.assertEquals(
				new  Biostar234081().instanceMain(new String[] {
					"-o",out.toString(),
					samFile
					}),0);
		
		support.assertIsValidBam(out);
		} finally {
			this.support.removeTmpFiles();
		}
	}
}
