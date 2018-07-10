package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar81455Test extends TestUtils {

@Test
public void test01() throws IOException {
	final File out = createTmpFile(".tsv");
	final File in = createTmpFile(".txt");
	final PrintWriter pw = new PrintWriter(in);
	pw.println("chr22\t41258261\nchr22\t52000000\nchr22\t0");
	pw.flush();
	pw.close();
	
	Assert.assertEquals(
		new Biostar81455().instanceMain(newCmd().
		add("-o").add(out).
		add("-KG").add(SRC_TEST_RESOURCE+"/test_vcf01.knownGenes.txt.gz").
		add(in).
		make()
		),0);
	
	super.assertTsvTableIsConsitent(out, null);
	
	}

}
