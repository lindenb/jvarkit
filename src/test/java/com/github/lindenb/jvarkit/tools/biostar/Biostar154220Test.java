package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.misc.SortSamRefName;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class Biostar154220Test {
	
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return (Object[][])support.
				allSamOrBams().
				map(F->new Object[] {F}).
				toArray()
				;
		}
	
	@Test(dataProvider="src1")
	public void test1(final String samFile) throws IOException {
		try {
		// limit to small ref
		if(SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(samFile)).
				getSequences().stream().
				mapToInt(L->L.getSequenceLength()).
				max().orElse(0) > 1_000_000) return;
		
		final Path out1 = support.createTmpPath(".bam");
		Assert.assertEquals(
				new  SortSamRefName().instanceMain(new String[] {
				"-o",out1.toString(),
				samFile
				}),0);
		
		
		support.assertIsValidBam(out1);
		
		final Path out2 = support.createTmpPath(".bam");
		Assert.assertEquals(
			new  Biostar154220().instanceMain(new String[] {
				"-o",out2.toString(),
				"-n","5",
				out1.toString()
			}),0);
		
		support.assertIsValidBam(out2);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}
}
