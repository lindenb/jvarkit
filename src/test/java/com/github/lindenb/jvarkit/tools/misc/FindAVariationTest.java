package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class FindAVariationTest extends TestUtils {

public void test01() throws IOException {
	final File input = createTmpFile(".tsv");
	PrintWriter pw=new PrintWriter(input);
	super._collectFiles(new File(SRC_TEST_RESOURCE),
			(D,F)->F.startsWith("S") && F.endsWith(".vcf.gz")
			).forEach(F->pw.println(F.getPath()));
	pw.flush();
	pw.close();
	Assert.assertTrue(wc(input)>0L);
	final File output = createTmpFile(".tsv");
	Assert.assertEquals(new FindAVariation().instanceMain(new String[]{
    		"-o",output.getPath(),
    		"-p","ref2:14",
    		input.getPath()
    		}),0);
	Assert.assertTrue(wc(output)>1L);
	super.assertTsvTableIsConsitent(output, null);
	}
}
