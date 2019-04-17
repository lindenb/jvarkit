package com.github.lindenb.jvarkit.util.picard;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class FourLinesFastqReaderTest
	{
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
		public Object[][] createData1() {
		 return new Object[][] {
			 {support.resource("SAMPLE1_GATGAATC_L002_R1_001.fastq.gz")},
			 {support.resource("S1.R1.fq.gz")},
		 };
		}

	
	 
	@Test(dataProvider="src1")
	public void test1(final String input) throws IOException
		{
		try {
		final FourLinesFastqReader r=new FourLinesFastqReader(new File(input));
		while(r.hasNext())
			{
			FastqRecord rec = r.next();
			Assert.assertNotNull(rec);
			}
		r.close();
		} finally {
			support.removeTmpFiles();
		}
		}

	@Test
	public void test2() throws Exception
		{
		try {
		final String fastqs=
				"@IL31_4368:1:1:996:8507/2\n" + 
				"\n" + 
				"+\n" + 
				"\n" + 
				"@IL31_4368:1:1:996:21421/2\n" + 
				"CAAAAACTTTCACTTTACCTGCCGGGTTTCCCAGTTTACATTCCACTGTTTGAC\n" + 
				"+\n" + 
				">DBDDB,B9BAA4AAB7BB?7BBB=91;+*@;5<87+*=/*@@?9=73=.7)7*\n" 
				;

		FourLinesFastqReader r=new FourLinesFastqReader(new ByteArrayInputStream(fastqs.getBytes()));
		r.setValidationStringency(ValidationStringency.LENIENT);
		Assert.assertTrue(r.hasNext());
		FastqRecord rec=r.next();
		Assert.assertNotNull(rec);
		Assert.assertTrue(rec.getBaseQualityString().isEmpty());
		Assert.assertTrue(rec.getReadString().isEmpty());
		Assert.assertTrue(r.hasNext());
		rec=r.next();
		Assert.assertNotNull(rec);
		Assert.assertFalse(rec.getBaseQualityString().isEmpty());
		Assert.assertFalse(rec.getReadString().isEmpty());
		Assert.assertFalse(r.hasNext());
		r.close();
		} finally {
			support.removeTmpFiles();
		}
		}
	}