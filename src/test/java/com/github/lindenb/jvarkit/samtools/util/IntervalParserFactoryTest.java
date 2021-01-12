package com.github.lindenb.jvarkit.samtools.util;

import java.nio.file.Paths;
import java.util.Optional;
import java.util.function.Function;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.StringUtilsTest;
import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtilsTest;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverterTest;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

@AlsoTest({SequenceDictionaryUtilsTest.class,SimpleIntervalTest.class,StringUtilsTest.class,ContigNameConverterTest.class})
public class IntervalParserFactoryTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void withoutDict() {
		Function<String,Optional<SimpleInterval>> parser = IntervalParserFactory.newInstance().make();
		final Optional<SimpleInterval> r=parser.apply("chr1:1-100");
		Assert.assertTrue(r.isPresent());
		Assert.assertEquals(r.get().getContig(),"chr1");
		Assert.assertEquals(r.get().getStart(),1);
		Assert.assertEquals(r.get().getEnd(),100);
		}

	@Test
	public void withComma() {
		Function<String,Optional<SimpleInterval>> parser = IntervalParserFactory.newInstance().make();
		final Optional<SimpleInterval> r=parser.apply("chr1:1-1,000");
		Assert.assertTrue(r.isPresent());
		Assert.assertEquals(r.get().getContig(),"chr1");
		Assert.assertEquals(r.get().getStart(),1);
		Assert.assertEquals(r.get().getEnd(),1000);
		}
	
	@Test
	public void withDict() {
		SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(support.resource("human_b37.dict")));
		Assert.assertNotNull(dict);
		Function<String,Optional<SimpleInterval>> parser = IntervalParserFactory.newInstance().
				dictionary(dict).make();
		final Optional<SimpleInterval> r=parser.apply("chr1:1-100");
		Assert.assertTrue(r.isPresent());
		Assert.assertEquals(r.get().getContig(),"1");
		Assert.assertEquals(r.get().getStart(),1);
		Assert.assertEquals(r.get().getEnd(),100);
		}
	
	@Test
	public void wholeContig() {
		SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(support.resource("human_b37.dict")));
		Assert.assertNotNull(dict);
		SAMSequenceRecord ssr=dict.getSequence(0);
		Function<String,Optional<SimpleInterval>> parser = IntervalParserFactory.newInstance().
				dictionary(dict).
				enableWholeContig().
				make();
		final Optional<SimpleInterval> r=parser.apply(ssr.getSequenceName());
		Assert.assertTrue(r.isPresent());
		Assert.assertEquals(r.get().getContig(),ssr.getSequenceName());
		Assert.assertEquals(r.get().getStart(),1);
		Assert.assertEquals(r.get().getEnd(),ssr.getSequenceLength());
		}
}
