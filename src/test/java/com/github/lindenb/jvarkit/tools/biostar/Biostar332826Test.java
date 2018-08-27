package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar332826Test extends TestUtils{
		
	@Test(dataProvider="all-vcf-files")
	public void test01(final String vcfpath) throws Exception
		{
		final File out = super.createTmpFile(".vcf");
		final File tmp = super.createTmpFile(".txt");
		final PrintWriter pw = new PrintWriter(tmp);
		for(int i=0;i< 1000;++i)
			{	
			pw.println("rs"+this.random.nextInt(10000));	
			}
		pw.flush();
		pw.close();
		
		Assert.assertEquals(new Biostar332826().instanceMain(new String[] {
			"-o",out.getPath(),
			"--ids",tmp.getPath(),
			vcfpath
			}),0);
		super.assertIsVcf(out);
		}
	
	}
