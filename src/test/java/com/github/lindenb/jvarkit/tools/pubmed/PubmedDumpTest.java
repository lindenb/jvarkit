package com.github.lindenb.jvarkit.tools.pubmed;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class PubmedDumpTest {
	
	protected final TestSupport support = new TestSupport();

	
@DataProvider(name = "src1")
public Object[][] getQueries() {
	return new Object[][] {
		{"Lindenbaum CloneIt"}
		};
	}

protected Path dumpAsXml(final String query ) throws IOException {
	final Path out = support.createTmpPath(".xml");
	Assert.assertEquals(
		new PubmedDump().instanceMain(new String[] {
		"-o",out.toString(),
		query
		}),0);
	support.assertIsXml(out);
	return out;
	}

@Test(dataProvider="src1")
public void test01(final String query ) throws IOException {
	final Path out = dumpAsXml(query);
	support.assertIsXml(out);
	}
}
