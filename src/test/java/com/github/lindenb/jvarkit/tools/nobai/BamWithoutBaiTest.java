package com.github.lindenb.jvarkit.tools.nobai;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.samtools.BamRecordGuesserTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.util.CoordMath;

@AlsoTest({LauncherTest.class,BamRecordGuesserTest.class})
public class BamWithoutBaiTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() 
		throws IOException
		{
		try {
			final Path out = support.createTmpPath(".bam");
			Assert.assertEquals(new BamWithoutBai().instanceMain(new String[] {
				"-r","chr3:38400000-38649667",
				"-o",out.toString(),
				"https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam"
				}),0);
			support.assertIsValidBam(out);
			
			Assert.assertTrue(support.samStream(out).allMatch(R->R.getReferenceName().equals("chr3") && CoordMath.overlaps(R.getStart(),R.getEnd(), 38400000,38649667)));
		} finally
			{
			support.removeTmpFiles();
			}
		}
}
