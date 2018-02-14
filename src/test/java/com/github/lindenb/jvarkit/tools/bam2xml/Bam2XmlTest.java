package com.github.lindenb.jvarkit.tools.bam2xml;

import java.io.File;
import java.io.IOException;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class Bam2XmlTest extends TestUtils {
@DataProvider(name = "src1")
public Object[][] createData1() {
	return new ParamCombiner().
		initList(collectAllSamOrBam()).
		build();
		}
	
@Test(dataProvider="src1")
public void test1(final String inBam) throws IOException {
	
	final File out = createTmpFile(".xml");
	final Bam2Xml cmd =new Bam2Xml();
	Assert.assertEquals(cmd.instanceMain(new String[] {
		"-o",out.getPath(),
		inBam
		}),0);
	assertIsXml(out);
	}
}
