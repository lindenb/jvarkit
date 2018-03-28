package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class SimplePlotTest extends TestUtils {
	
private File histo1(final boolean sort_uniq) throws IOException {
	final File in = super.createTmpFile(".txt");
	final PrintWriter pw = new PrintWriter(in);
	variantStream(new File(SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz")).
	map(V->V.getContig()).
	collect(Collectors.groupingBy(Function.identity(),Collectors.counting())).
	forEach((K,V)->{
			pw.println(sort_uniq?
					String.format(" %06d %s", V,K):
					K+"\t"+V
					);
			}
			);
	pw.flush();
	pw.close();
	return in;
	}



@Test
public void testHistogram01() throws IOException{
	final File in = histo1(false);
	final File imgOut = super.createTmpFile(".png");
	SimplePlot.main(new String[] {
			"-t","SIMPLE_HISTOGRAM",
			"-o",imgOut.getPath(),
			in.getPath()
			});
	super.assertIsImage(imgOut);
	}
@Test
public void testHistogram02() throws IOException{
	final File in = histo1(true);
	final File imgOut = super.createTmpFile(".png");
	SimplePlot.main(new String[] {
			"-t","SIMPLE_HISTOGRAM",
			"-su",
			"-o",imgOut.getPath(),
			in.getPath()
			});
	super.assertIsImage(imgOut);
	}
@Test
public void testPie01() throws IOException{
	final File in = histo1(false);
	final File imgOut = super.createTmpFile(".png");
	SimplePlot.main(new String[] {
			"-t","PIE",
			"-o",imgOut.getPath(),
			in.getPath()
			});
	super.assertIsImage(imgOut);
	}
@Test
public void testPie02() throws IOException{
	final File in = histo1(true);
	final File imgOut = super.createTmpFile(".png");
	SimplePlot.main(new String[] {
			"-t","PIE",
			"-su",
			"-o",imgOut.getPath(),
			in.getPath()
			});
	super.assertIsImage(imgOut);
	}

}
