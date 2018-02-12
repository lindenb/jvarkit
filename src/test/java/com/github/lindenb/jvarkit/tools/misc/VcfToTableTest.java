package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class VcfToTableTest  extends TestUtils {

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectAllVcfs()).
			build();
	}

@Test(dataProvider="src1")
public void test01(final String inputFile) 
	throws IOException
	{
	final File out = File.createTempFile(".tmp.", ".txt");
	final VcfToTable cmd =new VcfToTable();
	Assert.assertEquals(0,cmd.instanceMain(new String[] {
		"-o",out.getPath(),
		inputFile
		}));
	Assert.assertTrue(out.delete());
	}

	
}
