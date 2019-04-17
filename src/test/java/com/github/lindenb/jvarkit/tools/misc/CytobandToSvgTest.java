package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class CytobandToSvgTest {
	private final TestSupport support = new TestSupport();

@Test
public void test01() throws IOException {
	try {
		final Path cytobandFile = support.createTmpPath(".txt");
		PrintWriter pw = new PrintWriter(Files.newBufferedWriter(cytobandFile));
		pw.println("chr1\t0\t100\tp01\tgneg");
		pw.println("chr1\t100\t200\tq01\tgneg");
		pw.println("chr2\t0\t100\tp01\tgneg");
		pw.println("chr2\t100\t180\tq01\tgneg");
		pw.flush();
		pw.close();
		support.assertTsvTableIsConsitent(cytobandFile, null);
		
		final Path out = support.createTmpPath(".svg");
	
		 Assert.assertEquals(new CytobandToSvg().instanceMain(new String[]{
	     		"-o",out.toString(),
	     		"-C",cytobandFile.toString()
	     	}),0);
		 support.assertIsXml(out);
		} finally {
		support.removeTmpFiles();
	}
	}
}
