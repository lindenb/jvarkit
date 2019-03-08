package com.github.lindenb.jvarkit.util.vcf;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.StringUtilsTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

@AlsoTest({StringUtilsTest.class})
public class AFExtractorFactoryTest {

private final TestSupport support = new TestSupport();

	
@DataProvider(name="data01")
public Object[][] data01(){
	return new Object[][] {
		{"AC/AN",1,1},
		{"AF",1,1},
		{"x_*",1,0},
		{"x_*;,,AC/AN;;AF,,,,",3,2},
		{"",0,0},
		{";,,,;",0,0},
		{"*;AC/AN",1,1}
		};
	}
	
@Test(dataProvider="data01")
void test01(String key,int expect,int expectValids) {
	final AFExtractorFactory ex = new AFExtractorFactory();
	Assert.assertEquals(ex.parseFieldExtractors(key).size(), expect);
	}

@Test(dataProvider="data01")
void testValid(String key,int expect,int expectValids) throws IOException
	{
	final VCFFileReader ft =  new VCFFileReader(Paths.get(support.resource("test_vcf01.vcf")),false);
	final AFExtractorFactory ex = new AFExtractorFactory();
	final VCFHeader header =  ft.getFileHeader();
	
	final List<AFExtractorFactory.AFExtractor> L=ex.parseFieldExtractors(key);
	int n_valid =  0;
	for(final AFExtractorFactory.AFExtractor x:L)
		{
		n_valid+= x.validateHeader(header)?1:0;
		}
	Assert.assertEquals(n_valid,expectValids);
	ft.iterator().stream().forEach(V->{
		for(final AFExtractorFactory.AFExtractor x:L)
			{
			final List<Double> ld = x.parse(V);
			Assert.assertNotNull(ld);
			Assert.assertEquals(ld.size(),V.getAlternateAlleles().size());
			}
	});
	
	ft.close();
	}
}
