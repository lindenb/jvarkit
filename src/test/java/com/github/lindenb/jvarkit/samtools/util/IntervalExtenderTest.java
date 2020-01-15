package com.github.lindenb.jvarkit.samtools.util;

import java.io.IOException;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class IntervalExtenderTest {
	private final TestSupport support = new TestSupport();

	private void run(SAMSequenceDictionary dict) {
		final String contig="RF01";
		
		IntervalExtender x = IntervalExtender.of(dict, "0");
		Assert.assertNotNull(x);
		Assert.assertFalse(x.isChanging());
		Assert.assertFalse(IntervalExtender.of(dict, "100%").isChanging());
		Assert.assertFalse(IntervalExtender.of(dict, "0,0,0kb").isChanging());
		Assert.assertFalse(IntervalExtender.of(dict, "1.0").isChanging());
		
		x = IntervalExtender.of(dict, "3bp");
		Assert.assertTrue(x.isChanging());
		Assert.assertFalse(x.isExtendingByFraction());
		SimpleInterval r1=new SimpleInterval(contig,10,20);
		SimpleInterval r2=x.apply(r1);
		Assert.assertEquals(r1.getContig(),r2.getContig());
		Assert.assertEquals(r2.getStart(),7);
		Assert.assertEquals(r2.getEnd(),23);

		x = IntervalExtender.of(dict, "1");
		Assert.assertFalse(x.isExtendingByFraction());
		Assert.assertTrue(x.isChanging());
		r1=new SimpleInterval(contig,10,20);
		r2=x.apply(r1);
		Assert.assertEquals(r1.getContig(),r2.getContig());
		Assert.assertEquals(r2.getStart(),9);
		Assert.assertEquals(r2.getEnd(),21);
		
		x = IntervalExtender.of(dict, "-1");
		Assert.assertFalse(x.isExtendingByFraction());
		Assert.assertTrue(x.isChanging());
		r1=new SimpleInterval(contig,10,20);
		r2=x.apply(r1);
		Assert.assertEquals(r1.getContig(),r2.getContig());
		Assert.assertEquals(r2.getStart(),11);
		Assert.assertEquals(r2.getEnd(),19);
		
		
		x = IntervalExtender.of(dict, "200%");
		Assert.assertTrue(x.isChanging());
		Assert.assertTrue(x.isExtendingByFraction());
		r1=new SimpleInterval(contig,10,19);
		Assert.assertEquals(r1.getLengthOnReference(),10);
		r2=x.apply(r1);
		Assert.assertEquals(r1.getContig(),r2.getContig());
		Assert.assertEquals(r2.getStart(),4);
		Assert.assertEquals(r2.getEnd(),23);
		Assert.assertEquals(r2.getLengthOnReference(),20);

		
		x = IntervalExtender.of(dict, "1");
		Assert.assertFalse(x.isExtendingByFraction());
		Assert.assertTrue(x.isChanging());
		r1=new SimpleInterval(contig,10,10);
		r2=x.apply(r1);
		Assert.assertEquals(r1.getContig(),r2.getContig());
		Assert.assertEquals(r2.getStart(),9);
		Assert.assertEquals(r2.getEnd(),11);

		x = IntervalExtender.of(dict, "150%");
		Assert.assertTrue(x.isExtendingByFraction());
		Assert.assertTrue(x.isChanging());
		r1=new SimpleInterval(contig,10,10);
		r2=x.apply(r1);
		Assert.assertEquals(r1.getContig(),r2.getContig());
		Assert.assertEquals(r2.getStart(),9);
		Assert.assertEquals(r2.getEnd(),10);

		
	}
	
	@Test
	public void test01() throws IOException {
		run(SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(support.resource("rotavirus_rf.dict"))));
	}
}
