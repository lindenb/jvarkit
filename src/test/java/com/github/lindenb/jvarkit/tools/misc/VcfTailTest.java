package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfTailTest {

	private final TestSupport support = new TestSupport();

	@DataProvider(name="src1")
	public Object[][] testData01() {
			return support.combine2(
					support.allVcfOrBcf(),
					Arrays.stream(new Object[] {0,1,5,1000})
					);
			}

@Test(dataProvider="src1")
public void test01(final String inputFile,int num) 
	throws IOException
	{
	try {
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(0,new VcfTail().instanceMain(new String[] {
			"-n",String.valueOf(num),
			"-o",out.toString(),
			inputFile
			}));
		Assert.assertTrue(support.variantStream(out).count() <=num);
		}
	finally {
		support.removeTmpFiles();
		}
	}

	
}
