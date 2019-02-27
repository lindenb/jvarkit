package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class Biostar3654Test {
	
private final TestSupport support = new TestSupport();
	
public void test01() throws IOException
	{
	try {
		final Path out = support.createTmpPath(".txt");
		Assert.assertEquals(
				new Biostar59647().instanceMain(new String[] {
				"-o",out.toString(),
				support.resource("rotavirus_rf.blastn.01.xml")
				}),0);
		support.assertIsNotEmpty(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}
}
