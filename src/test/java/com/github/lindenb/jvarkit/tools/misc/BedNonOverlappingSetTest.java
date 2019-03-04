package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


import org.testng.Assert;
import org.testng.annotations.Test;

public class BedNonOverlappingSetTest {
	
	private final TestSupport support = new TestSupport();

	
	@Test
    public void test01() throws IOException {
	try {
	final Path tmpBed =support.createTmpPath(".bed");
    final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(tmpBed));
    pw.println("RF01\t1\t100");
    pw.println("RF01\t10\t110");

	pw.flush();
	pw.close();
	support.assertIsBed(tmpBed);
	Path pat= Files.createTempFile("tmp.__SETID__.", ".bed");
	Assert.assertEquals(
			new BedNonOverlappingSet().instanceMain(new String[] {
					"-o",pat.toString(),
					"-x","10",
					"-R",support.resource("rotavirus_rf.fa"),
					tmpBed.toString()
			}),0);
	for(final Path bf: Files.list(pat.getParent()).filter(F->F.getFileName().toString().startsWith("tmp.") && F.getFileName().toString().endsWith(".bed")).collect(Collectors.toList()))
		{
		support.assertIsBed(bf);
		Files.deleteIfExists(bf);
		}
	
	Files.deleteIfExists(pat);
	} finally {
		support.removeTmpFiles();
	}
	}
}
