package com.github.lindenb.jvarkit.util.samtools;

import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class ContigDictComparatorTest
	{
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
	 return new Object[][] {
		 {support.resource("rotavirus_rf.dict")},
		 {support.resource("toy.dict")},
		 {support.resource("human_b37.dict")}
	 };
	}
	
	@Test(dataProvider="src1")
	public void test(final String dictfile) {
		try {
		SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(dictfile));
		Assert.assertNotNull(dict);
		List<String> contigs = new ArrayList<>();
		while(contigs.size()<1000) {
			contigs.add(dict.getSequence(support.random.nextInt(dict.size())).getSequenceName());
		}
		Collections.shuffle(contigs, support.random);
		Collections.sort(contigs, new ContigDictComparator(dict));
		
		for(int i=1;i< contigs.size();i++)
			{
			Assert.assertTrue(dict.getSequenceIndex(contigs.get(i-1))<=dict.getSequenceIndex(contigs.get(i)));
			}
		} finally {
			support.removeTmpFiles();
		}
		}
	}