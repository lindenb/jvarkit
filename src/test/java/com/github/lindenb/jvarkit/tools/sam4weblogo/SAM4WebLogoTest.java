package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;


public class SAM4WebLogoTest extends TestUtils
	{
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectIndexedBams()).
			build();
		}
	@Test(dataProvider="src1")
	public void test1(final String inBam) throws IOException {
		final File out = createTmpFile(".txt");
		final SAMSequenceDictionary dict= SAMSequenceDictionaryExtractor.extractDictionary(new File(inBam));
		if(dict==null || dict.isEmpty()) return;
		final SAMSequenceRecord ssr = dict.getSequence(this.random.nextInt(dict.size()));
		final int pos= this.random.nextInt(ssr.getSequenceLength());
		Assert.assertEquals(0,new SAM4WebLogo().instanceMain(newCmd().add(
        		"-o",out.getPath(),
        		"-r",ssr.getSequenceName()+":"+pos+"-"+Math.min(pos+this.random.nextInt(ssr.getSequenceLength()),ssr.getSequenceLength()),
        		inBam
        		).addIf(random.nextBoolean(),"-c").
				make()));
		assertIsNotEmpty(out);
		}
	}
