package com.github.lindenb.jvarkit.util.vcf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.IOUtilsTest;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.CharSplitterTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.tools.vcftrios.DeNovoDetectorTest;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtilsTest;

@AlsoTest({IOUtilsTest.class,CharSplitterTest.class,
	SequenceDictionaryUtilsTest.class,
	AFExtractorFactoryTest.class,
	DeNovoDetectorTest.class
	})
public class VCFUtilsTest {
final private TestSupport support= new TestSupport();
@DataProvider(name="src01")
public Object[][] testData01() {
	return new Object[][] {
		{"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"}
	};
}

	
@Test(dataProvider="src01")
public void test01(String header) throws IOException {
	VCFUtils.parseHeader(new BufferedReader(new StringReader(header)));
	}
@Test(dataProvider="src01")
public void test02(String header) throws IOException {
	VCFUtils.parseHeader(Arrays.asList(CharSplitter.NEWLINE.split(header)));
	}
@Test(dataProvider="src01")
public void test03(String header) throws IOException {
	VCFUtils.parseHeader(IOUtils.toLineIterator(new BufferedReader(new StringReader(header))));
	}
@Test
public void testIsAVcf() {
	Path p = Paths.get(support.resource("rotavirus_rf.ann.vcf.gz"));
	Assert.assertTrue(VCFUtils.isVcfPath(p));
	Assert.assertTrue(VCFUtils.isTabixVcfPath(p));
	Assert.assertFalse(VCFUtils.isTribbleVcfPath(p));
	
	File f = new File(support.resource("rotavirus_rf.ann.vcf.gz"));
	Assert.assertTrue(VCFUtils.isVcfFile(f));
	Assert.assertTrue(VCFUtils.isTabixVcfFile(f));
	Assert.assertFalse(VCFUtils.isTribbleVcfFile(f));
	
	p = Paths.get(support.resource("test_vcf01.vcf"));
	Assert.assertTrue(VCFUtils.isVcfPath(p));
	Assert.assertFalse(VCFUtils.isTabixVcfPath(p));
	Assert.assertTrue(VCFUtils.isTribbleVcfPath(p));

	f = new File(support.resource("test_vcf01.vcf"));
	Assert.assertTrue(VCFUtils.isVcfFile(f));
	Assert.assertFalse(VCFUtils.isTabixVcfFile(f));
	Assert.assertTrue(VCFUtils.isTribbleVcfFile(f));

	
	p = Paths.get(support.resource("toy.bam"));
	Assert.assertFalse(VCFUtils.isVcfPath(p));
	Assert.assertFalse(VCFUtils.isTabixVcfPath(p));
	Assert.assertFalse(VCFUtils.isTribbleVcfPath(p));

	}

}
