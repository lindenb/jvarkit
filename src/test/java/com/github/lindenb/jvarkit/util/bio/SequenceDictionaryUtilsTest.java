package com.github.lindenb.jvarkit.util.bio;

import java.io.File;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class SequenceDictionaryUtilsTest  extends TestUtils  {
@Test
public void testIsGRCh37() {
	Assert.assertTrue(SequenceDictionaryUtils.isGRCh37(
			SequenceDictionaryUtils.extractRequired(
					new File(SRC_TEST_RESOURCE+"/human_b37.dict"))));
	
	Assert.assertFalse(SequenceDictionaryUtils.isGRCh37(
			SequenceDictionaryUtils.extractRequired(
					new File(SRC_TEST_RESOURCE+"/rotavirus_rf.fa"))));
	}
@Test
public void testIsHuman() {
	Assert.assertTrue(SequenceDictionaryUtils.isHuman(
			SequenceDictionaryUtils.extractRequired(
					new File(SRC_TEST_RESOURCE+"/human_b37.dict"))));
	
	Assert.assertFalse(SequenceDictionaryUtils.isHuman(
			SequenceDictionaryUtils.extractRequired(
					new File(SRC_TEST_RESOURCE+"/rotavirus_rf.fa"))));
	}

}
