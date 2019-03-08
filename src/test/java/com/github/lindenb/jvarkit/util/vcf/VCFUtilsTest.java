package com.github.lindenb.jvarkit.util.vcf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.IOUtilsTest;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.CharSplitterTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtilsTest;

@AlsoTest({IOUtilsTest.class,CharSplitterTest.class,
	SequenceDictionaryUtilsTest.class,
	AFExtractorFactoryTest.class
	})
public class VCFUtilsTest {

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

}
