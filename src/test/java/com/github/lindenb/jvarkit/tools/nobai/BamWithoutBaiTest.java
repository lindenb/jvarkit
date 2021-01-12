package com.github.lindenb.jvarkit.tools.nobai;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.samtools.BamRecordGuesserTest;
import com.github.lindenb.jvarkit.tests.AlsoTest;
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
		final String contig= "chr3";
		final int chromStart=38548061;
		final int chromEnd=38649667;
		try {
			// tmp output bam
			final Path out = support.createTmpPath(".bam");
			// execute
			Assert.assertEquals(new BamWithoutBai().instanceMain(new String[] {
				"-r",contig+":"+chromStart+"-"+chromEnd,
				"-o",out.toString(),
				"https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam"
				}),0);
			// check this is bam
			support.assertIsValidBam(out);
			//check overlap
			Assert.assertTrue(support.samStream(out).allMatch(R->R.getReferenceName().equals(contig) && CoordMath.overlaps(R.getStart(),R.getEnd(), chromStart,chromEnd)));
			// check first read
			Assert.assertEquals(support.samStream(out).map(R->R.getReadName()).findFirst().orElse(null),"D2FC08P1:268:C3NPCACXX:8:1111:11204:70951");
			// check last read
			Assert.assertEquals(support.samStream(out).map(R->R.getReadName()).reduce((A,B)->B).orElse(null),"D2FC08P1:268:C3NPCACXX:8:2312:16447:12679");
		} finally
			{
			support.removeTmpFiles();
			}
		}
}
