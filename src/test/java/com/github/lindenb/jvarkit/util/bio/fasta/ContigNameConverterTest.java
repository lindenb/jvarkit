package com.github.lindenb.jvarkit.util.bio.fasta;

import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.StringUtilsTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtilsTest;

import htsjdk.samtools.SAMSequenceDictionary;

@AlsoTest({StringUtilsTest.class,SequenceDictionaryUtilsTest.class})
public class ContigNameConverterTest {
	private  final TestSupport support = new TestSupport();

	@Test
	public void testFromOneDictionary() {
		final SAMSequenceDictionary dict =SequenceDictionaryUtils.extractRequired(Paths.get(support.resource("human_b37.dict")));
		final ContigNameConverter f= ContigNameConverter.fromOneDictionary(dict);
		Assert.assertEquals(f.apply("1"), "1");
		Assert.assertNull(f.apply("z"));
		Assert.assertNotEquals(f.apply("11"), "1");
		Assert.assertEquals(f.apply("chr1"), "1");
		Assert.assertNotEquals(f.apply("chr1"), "chr1");
		Assert.assertEquals(f.apply("chrX"), "X");
		Assert.assertEquals(f.apply("chrY"), "Y");
		Assert.assertEquals(f.apply("chrM"), "MT");
		}
}
