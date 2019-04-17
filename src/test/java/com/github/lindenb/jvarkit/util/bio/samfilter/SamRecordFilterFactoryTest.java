package com.github.lindenb.jvarkit.util.bio.samfilter;

import java.io.IOException;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;

public class SamRecordFilterFactoryTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void testFilters() throws IOException {
		final SamRecordFilterFactory srff = new SamRecordFilterFactory();
		SamRecordFilter f1 = srff.buildDefault();
		SamRecordFilter f2 = srff.convert("");
		SamRecordFilter f3 = srff.convert("isMapped()");
		
		SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(Paths.get(support.resource("S1.bam")));
		SAMRecordIterator iter = sr.iterator();
		while(iter.hasNext()) {
			final SAMRecord rec=iter.next();
			Assert.assertFalse(f2.filterOut(rec));
			if(rec.getReadUnmappedFlag()) Assert.assertTrue(f1.filterOut(rec));
			if(!rec.getReadUnmappedFlag()) Assert.assertTrue(f3.filterOut(rec));
		}
		iter.close();
		sr.close();
	}
}