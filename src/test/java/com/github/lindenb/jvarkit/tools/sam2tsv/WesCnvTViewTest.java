package com.github.lindenb.jvarkit.tools.sam2tsv;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class WesCnvTViewTest  extends TestUtils {
	@Test
	public void test01() throws IOException {
		final  File bamList = super.createTmpFile(".list");
		final PrintWriter pw = new PrintWriter(bamList);
		for(int i=1;i<=5;i++)
			{
			pw.println(SRC_TEST_RESOURCE+"/S"+i+".bam");
			}
		pw.flush();
		pw.close();
		final File out = super.createTmpFile(".txt");
		Assert.assertEquals(new WesCnvTView().instanceMain(new String[] {
				"-l",bamList.getPath(),
				"-o",out.getPath(),
				"RF01:200-1000",
				"RF02:300-500"
				}),0);
		assertIsNotEmpty(out);
	}
}
