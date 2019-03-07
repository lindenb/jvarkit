package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.testng.annotations.Test;
import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class SimplePlotTest {
	private final TestSupport support = new TestSupport();
	
	
	private Path histo1(final boolean sort_uniq) throws IOException {
	final Path in = support.createTmpPath(".txt");
	final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(in));
	support.variantStream(Paths.get(support.resource("rotavirus_rf.vcf.gz"))).
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



@Test(enabled=false)
public void testHistogram01() throws IOException{
	try {
	final Path in = histo1(false);
	final Path imgOut = support.createTmpPath(".R");
	launch(new String[] {
			"-t","SIMPLE_HISTOGRAM",
			"-o",imgOut.toString(),
			in.toString()
			});
} finally {
	support.removeTmpFiles();
}
	}
@Test(enabled=false)
public void testHistogram02() throws IOException{
	try {
	final Path in = histo1(true);
	final Path imgOut = support.createTmpPath(".R");
	launch(new String[] {
			"-t","SIMPLE_HISTOGRAM",
			"-su",
			"-o",imgOut.toString(),
			in.toString()
			});
} finally {
	support.removeTmpFiles();
}
	}
@Test(enabled=false)
public void testPie01() throws IOException{
	try {
	final Path in = histo1(false);
	final Path imgOut = support.createTmpPath(".R");
	launch(new String[] {
			"-t","PIE",
			"-o",imgOut.toString(),
			in.toString()
			});
		} finally {
			support.removeTmpFiles();
		}
	}
@Test(enabled=false)
public void testPie02() throws IOException{
	try {
	final Path in = histo1(true);
	final Path imgOut = support.createTmpPath(".R");
	launch(new String[] {
			"-t","PIE",
			"-su",
			"-o",imgOut.toString(),
			in.toString()
			});
	} finally {
		support.removeTmpFiles();
	}
	}

private void launch(final String[] args) {
	Assert.assertEquals(new SimplePlot().instanceMain(args),0);
	}
}
