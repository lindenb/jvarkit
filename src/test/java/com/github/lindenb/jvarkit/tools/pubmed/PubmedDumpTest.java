package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class PubmedDumpTest extends TestUtils{
@DataProvider(name = "src1")
public Object[][] getQueries() {
	return new Object[][] {
		{"Lindenbaum CloneIt"}
		};
	}

@Test(dataProvider="src1")
public void test01(final String query ) throws IOException {
	final File out = createTmpFile(".xml");
	Assert.assertEquals(
		new PubmedDump().instanceMain(newCmd().
		add("-o").add(out).
		add(query).
		make()
		),0);
	assertIsXml(out);
	}
}
