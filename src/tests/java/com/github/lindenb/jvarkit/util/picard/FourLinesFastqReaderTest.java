package com.github.lindenb.jvarkit.util.picard;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.StringReader;

import org.testng.Assert;
import org.testng.annotations.Test;


public class FourLinesFastqReaderTest
	{
	@Test()
	public void test1() throws Exception
		{
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
		
	
		}
	}