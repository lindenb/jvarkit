package com.github.lindenb.jvarkit.util.bio.bed;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.CharSplitterTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;

@AlsoTest({CharSplitterTest.class})
public class BedLineCodecTest {
	@Test
	public void test() {
		BedLineCodec codec = new BedLineCodec();
		Assert.assertNull(codec.decode(""));
		Assert.assertNull(codec.decode("#"));
		Assert.assertNull(codec.decode("browser"));
		Assert.assertNull(codec.decode("track"));
		final BedLine bed =codec.decode("chr1\t1\t100");
		Assert.assertNotNull(bed);
		Assert.assertEquals(bed.getContig(), "chr1");
		Assert.assertEquals(bed.getStart(),2);
		Assert.assertEquals(bed.getEnd(),100);
		Assert.assertEquals(bed.get(1),"1");
		Assert.assertEquals(bed.getOrDefault(2,"z"),"100");
		Assert.assertEquals(bed.getOrDefault(3,"z"),"z");
	}
}
