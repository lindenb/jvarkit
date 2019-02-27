package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Biostar145820Test{
	
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
	public void testShuffle(final String bam) throws IOException {
		try {
			final Path out = support.createTmpPath(".bam");
			Assert.assertEquals(
				new Biostar145820().instanceMain(new String[] {
				"-o",out.toString(),
				bam})
				,0);
			support.assertIsValidBam(out);
			Assert.assertEquals(support.wc(out), support.wc(Paths.get(bam)));
			} 
		finally
			{
			support.removeTmpFiles();
			}
		}
	
	@Test(dataProvider="src1")
	public void testDownSample(final String bam) throws IOException {
		try {
			final Path out = support.createTmpPath(".bam");
			Assert.assertEquals(
				new Biostar145820().instanceMain(new String[] {
				"-o",out.toString(),
				"-n","10",
				bam})
				,0);
			support.assertIsValidBam(out);
			Assert.assertTrue(support.wc(out)<=10);
			} 
		finally
			{
			support.removeTmpFiles();
			}
		}
	
	}
