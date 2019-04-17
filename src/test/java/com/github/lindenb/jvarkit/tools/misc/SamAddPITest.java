package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

public class SamAddPITest {
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
		// limit to small ref
		if(SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(samFile)).
				getSequences().stream().
				mapToInt(L->L.getSequenceLength()).
				min().orElse(0) > 1_000_000) return;

		
		final Path out = support.createTmpPath(".bam");
		final Path in = support.addClippingToBam(Paths.get(samFile));
		Assert.assertEquals(
			new SamAddPI().instanceMain(new String[] {
					"-o",out.toString(),"-w","-N","1000",
					in.toString()
			}),0);
		support.assertIsValidBam(out);
		}
		finally {
			support.removeTmpFiles();
		}
	}
}
