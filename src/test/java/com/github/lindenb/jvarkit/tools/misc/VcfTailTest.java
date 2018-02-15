package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfTailTest  extends TestUtils {

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectAllVcfs()).
			product(0,1,5,1000).
			build();
	}

@Test(dataProvider="src1")
public void test01(final String inputFile,int num) 
	throws IOException
	{
	final File out = super.createTmpFile(".vcf");
	final VcfTail cmd =new VcfTail();
	Assert.assertEquals(0,cmd.instanceMain(new String[] {
		"-n",String.valueOf(num),
		"-o",out.getPath(),
		inputFile
		}));
	Assert.assertTrue(variantStream(out).count() <=num);
	
	}

	
}
