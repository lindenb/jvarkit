package com.github.lindenb.jvarkit.samtools;

import java.io.File;
import java.io.IOException;
import java.io.PushbackInputStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloserUtil;

public class BamRecordGuesserTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void testFindSam() 
		throws IOException
		{
		BlockCompressedInputStream bcis = null;
		try {
			File samFile=new File(this.support.resource("S1.bam"));
			Assert.assertTrue(samFile.exists());
			bcis = new BlockCompressedInputStream(samFile);
			PushbackInputStream pbis = new PushbackInputStream(bcis,BamRecordGuesser.BUFFER_SIZE);
			BamRecordGuesser g  = new BamRecordGuesser(SamReaderFactory.make().getFileHeader(samFile));
			boolean b=g.find(pbis);
			bcis.close();
			bcis=null;
			Assert.assertTrue(b);
		} finally
			{
			CloserUtil.close(bcis);
			support.removeTmpFiles();
			}
		}

}
