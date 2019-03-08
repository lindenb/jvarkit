package com.github.lindenb.jvarkit.util.bio;

import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.StringUtilsTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverterTest;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

@AlsoTest({SequenceDictionaryUtilsTest.class,StringUtilsTest.class,ContigNameConverterTest.class})
public class IntervalParserTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void withDict() {
		SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(support.resource("human_b37.dict")));
		Assert.assertNotNull(dict);
		IntervalParser parser = new IntervalParser(dict);
		Assert.assertEquals(parser.parse("chr1:1-100"),new Interval("1",1,100));
		}
}
