package com.github.lindenb.jvarkit.util.samtools;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class ContigDictComparatorTest extends TestUtils
	{
	@Test(dataProvider="all-dictionaries")
	public void test(final String dictfile) {
		SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(new File(dictfile));
		Assert.assertNotNull(dict);
		List<String> contigs = new ArrayList<>();
		while(contigs.size()<1000) {
			contigs.add(dict.getSequence(super.random.nextInt(dict.size())).getSequenceName());
		}
		Collections.shuffle(contigs, this.random);
		Collections.sort(contigs, new ContigDictComparator(dict));
		
		for(int i=1;i< contigs.size();i++)
			{
			Assert.assertTrue(dict.getSequenceIndex(contigs.get(i-1))<=dict.getSequenceIndex(contigs.get(i)));
			}
		}
	}