package com.github.lindenb.jvarkit.tools.cmpbams;

import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class CommBamsTest{
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
	public void test1(final String bampath) throws Exception
		{
		try {
		final Path out = support.createTmpPath(".txt");
		final Path sortedBam1= support.sortBamOnQueryName(Paths.get(bampath),R->support.random.nextDouble()<0.5);
		final Path sortedBam2= support.sortBamOnQueryName(Paths.get(bampath),R->support.random.nextDouble()<0.5);
		Assert.assertEquals(new CommBams().instanceMain(new String[] {
			"-o",out.toString(),
			sortedBam1.toString(),
			sortedBam2.toString(),
			}),0);
		} finally {
			support.removeTmpFiles();
		}
		}
}
