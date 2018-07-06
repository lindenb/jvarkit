package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar3654Test extends TestUtils {
public void test01() throws IOException
	{
	final File out = super.createTmpFile(".txt");
	Assert.assertEquals(
			new Biostar59647().instanceMain(newCmd().
			add("-o").add(out).
			add(SRC_TEST_RESOURCE+"/rotavirus_rf.blastn.01.xml").
			make()
			),0);
	assertIsNotEmpty(out);
	}
}
