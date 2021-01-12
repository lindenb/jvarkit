package com.github.lindenb.jvarkit.tools.bioalcidae;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class BioAlcidaeJdkTest
	{
	private final TestSupport support = new TestSupport();

	
	

	@Test
	public void testVcf() 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[] {
				"-e","println(stream().count());",
				"-o",out.toString(),
				support.resource("toy.vcf.gz")
				}),0);
			support.assertIsNotEmpty(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		
		}

	
	

	@Test
	public void testBam() 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[] {
				"-e","println(stream().count());",
				"-o",out.toString(),
				support.resource("toy.bam")
				}),0);
			support.assertIsNotEmpty(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		
		}

	@Test
	public void testFasta() 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[] {
				"-e","println(stream().count());",
				"-o",out.toString(),
				support.resource("toy.fa")
				}),0);
			support.assertIsNotEmpty(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		
		}
	
	
	@Test
	public void testCram()  throws IOException
		{
		try {
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[] {
				"-e","println(stream().count());",
				"-o",out.toString(),
				"-R",support.resource("toy.fa"),
				support.resource("toy.cram")
				}),0);
			support.assertIsNotEmpty(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		
		}
	@Test
	public void testGtf()  throws IOException
		{
		try {
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[] {
				"-e","println(stream().count());",
				"-o",out.toString(),
				support.resource("Homo_sapiens.GRCh37.87.gtf.gz")
				}),0);
			support.assertIsNotEmpty(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
		
		}
}
