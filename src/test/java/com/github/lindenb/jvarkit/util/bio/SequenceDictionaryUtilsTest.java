package com.github.lindenb.jvarkit.util.bio;

import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.StringUtilsTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverterTest;

@AlsoTest({StringUtilsTest.class,ContigNameConverterTest.class,IntervalParserTest.class})
public class SequenceDictionaryUtilsTest  {
	private  final TestSupport support = new TestSupport();

@Test
public void testIsGRCh37() {
	Assert.assertTrue(SequenceDictionaryUtils.isGRCh37(
			SequenceDictionaryUtils.extractRequired(
					Paths.get(support.resource("human_b37.dict")))));
	
	Assert.assertFalse(SequenceDictionaryUtils.isGRCh37(
			SequenceDictionaryUtils.extractRequired(
					Paths.get(support.resource("rotavirus_rf.fa")))));
	
	Assert.assertFalse(SequenceDictionaryUtils.getBuildName(
			SequenceDictionaryUtils.extractRequired(
					Paths.get(support.resource("rotavirus_rf.fa")))).isPresent());
	
	Assert.assertTrue(SequenceDictionaryUtils.getBuildName(
			SequenceDictionaryUtils.extractRequired(
					Paths.get(support.resource("human_b37.dict")))).isPresent());

	}
@Test
public void testIsHuman() {
	Assert.assertTrue(SequenceDictionaryUtils.isHuman(
			SequenceDictionaryUtils.extractRequired(
					Paths.get(support.resource("human_b37.dict")))));
	
	Assert.assertFalse(SequenceDictionaryUtils.isHuman(
			SequenceDictionaryUtils.extractRequired(
					Paths.get(support.resource("rotavirus_rf.fa")))));
	}

}
