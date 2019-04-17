package com.github.lindenb.jvarkit.tools.biostar;

import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class Biostar332826Test {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return (Object[][])support.
				allVcfOrBcf().
				map(F->new Object[] {F}).
				toArray()
				;
		} 
	
	@Test(dataProvider="src1")
	public void test01(final String vcfpath) throws Exception
		{
		try {
			final Path out = support.createTmpPath(".vcf");
			final Path tmp = support.createTmpPath(".txt");
			final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(tmp));
			for(int i=0;i< 1000;++i)
				{	
				pw.println("rs"+support.random.nextInt(10000));	
				}
			pw.flush();
			pw.close();
			
			Assert.assertEquals(new Biostar332826().instanceMain(new String[] {
				"-o",out.toString(),
				"--ids",tmp.toString(),
				vcfpath
				}),0);
			support.assertIsVcf(out);
		} finally
			{
			support.removeTmpFiles();
			}
		}
	
	}
