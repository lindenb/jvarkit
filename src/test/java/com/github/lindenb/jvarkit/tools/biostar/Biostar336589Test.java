package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar336589Test extends TestUtils {
@Test
public void test01() throws IOException {
	final File bed = super.createTmpFile(".bed");
	final PrintWriter pw = new PrintWriter(bed);
	pw.println("RF01\t101\t1000\tf1\t100\t+\t.\t.\t0,255,0");
	pw.println("RF02\t101\t1000\tf1");
	pw.println("RF03\t101\t1000\tf1\t1000");
	pw.close();
	assertIsBed(bed);
	
	
	final File out = super.createTmpFile(".svg");
	Assert.assertEquals(new Biostar336589().instanceMain(new String[] {
			"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.dict",
			"-o",out.getPath(),
			bed.getPath(),
			}),0);
	assertIsXml(out);
	}
}
