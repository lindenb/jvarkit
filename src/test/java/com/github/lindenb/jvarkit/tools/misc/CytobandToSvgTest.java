package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class CytobandToSvgTest extends TestUtils{
@Test
public void test01() throws IOException {
	final File cytobandFile = super.createTmpFile(".txt");
	PrintWriter pw = new PrintWriter(cytobandFile);
	pw.println("chr1\t0\t100\tp01\tgneg");
	pw.println("chr1\t100\t200\tq01\tgneg");
	pw.println("chr2\t0\t100\tp01\tgneg");
	pw.println("chr2\t100\t180\tq01\tgneg");
	pw.flush();
	pw.close();
	assertTsvTableIsConsitent(cytobandFile, null);
	
	final File out = super.createTmpFile(".svg");

	 Assert.assertEquals(new CytobandToSvg().instanceMain(new String[]{
     		"-o",out.getPath(),
     		"-C",cytobandFile.getPath()
     	}),0);
	assertIsXml(out);
	}
}
