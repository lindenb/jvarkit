package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class ValidateCnvTest  extends TestUtils{
	@Test
	public void test01() throws IOException
		{
		final File bed= super.createTmpFile(".ped");
		final PrintWriter pw = new PrintWriter(bed);
		for(int i=0;i< 10;++i) {
			int p = this.random.nextInt(100);
			int L = this.random.nextInt(100);
			pw.println("RF01"+i+"\t"+p+"\t"+(p+L)+"\tid"+i);
			}
		pw.flush();
		pw.close();
		
		
		final File out = super.createTmpFile(".vcf");
		Assert.assertEquals(new ValidateCnv().instanceMain(new String[] {
				"-o",out.getPath(),
				"-B",bed.getPath(),
			SRC_TEST_RESOURCE+"/S1.bam",
			SRC_TEST_RESOURCE+"/S2.bam",
			SRC_TEST_RESOURCE+"/S3.bam",
			SRC_TEST_RESOURCE+"/S4.bam"
			}),0
			);
		assertIsVcf(out);
		}
	}
