package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class PubmedDumpTest {
	
	
@DataProvider(name = "src1")
public Object[][] getQueries() {
	return new Object[][] {
		{"Lindenbaum CloneIt"}
		};
	}

protected Path dumpAsXml(final TestSupport support ,final String query ) throws IOException {
	final Path out = support.createTmpPath(".xml");
	Assert.assertEquals(
		new PubmedDump().instanceMain(new String[] {
		"-o",out.toString(),
		query
		}),0,"query was "+query);
	support.assertIsXml(out);
	return out;
	}

@Test(dataProvider="src1")
public void test01(final String query ) throws IOException {
	final TestSupport support = new TestSupport();
	try {
		final Path out = dumpAsXml(support,query);
		support.assertIsXml(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}
}
