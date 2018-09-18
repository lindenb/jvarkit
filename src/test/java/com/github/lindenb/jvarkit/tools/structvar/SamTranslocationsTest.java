package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class SamTranslocationsTest  extends TestUtils {
	@Test(dataProvider="all-sam-or-bam-files")
	public void test01(final String inBam) 
		throws IOException
		{
		final File out = super.createTmpFile(".txt");
		Assert.assertEquals(new SamTranslocationsTest().instanceMain(new String[] {
			"-o",out.getPath(),
			inBam
			}),0);
		assertTsvTableIsConsitent(out, null);
		}
}
