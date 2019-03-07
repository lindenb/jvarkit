package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class SamToPslTest {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allSamOrBams().
				map(F->new Object[] {F})
				)
				;
		}
	
	@Test(dataProvider="src1")
	public void test1(final String samFile) throws IOException {
		try {
			final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(samFile));
			if(dict==null || dict.isEmpty()) return;
			if(dict.getSequences().stream().
					mapToInt(L->L.getSequenceLength()).
					max().orElse(0) > 1_000_000) return;
	 
			final Path out = support.createTmpPath(".txt");
			Assert.assertEquals(
				new SamToPsl().instanceMain(new String[] {
						"-o",out.toString(),
						samFile
				}),0);
			support.assertIsNotEmpty(out);
			} 
		finally
			{
			support.removeTmpFiles();
			}
	}
}
