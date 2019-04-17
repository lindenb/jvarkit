package com.github.lindenb.jvarkit.tools.biostar;

import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Biostar352930Test {
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
	public void test01(final String bampath) throws Exception
		{
		try {
			final Path out = support.createTmpPath(".bam");
			final Path sortedBam1= support.sortBamOnQueryName(Paths.get(bampath),null);
			
			Assert.assertEquals(new Biostar352930().instanceMain(new String[] {
				"-o",out.toString(),
				sortedBam1.toString()
				}),0);
			support.assertIsValidBam(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	}
