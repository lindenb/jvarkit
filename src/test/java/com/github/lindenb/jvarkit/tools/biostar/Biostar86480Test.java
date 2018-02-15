package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class Biostar86480Test extends TestUtils{
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new ParamCombiner().
			initList(collectAllFasta()).
			product("-E HhaI","-E AluI","-E HhaI -E AluI").
			build();
		}
		
	@Test(dataProvider="src1")
	public void test1(final String infFasta,final String enz) throws IOException {
		final File out = createTmpFile(".txt");
		final Biostar86480 cmd =new Biostar86480();
		Assert.assertEquals(
			cmd.instanceMain(newCmd().
			add("-o",out.getPath()).
			split(enz).
			add(infFasta).
			make()
			),0);
		assertIsNotEmpty(out);
		}
}
