package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class NcbiGeneDumpTest extends TestUtils{
@Test
public void test01() throws IOException {
	final File out =super.createTmpFile(".xml");
	Assert.assertEquals(new NcbiGeneDump().instanceMain(new String[] {
			"-o",out.getPath(),
			"SCN5A","NOTCH2"
			}),0);
	super.assertIsXml(out);
	}

@Test
public void testXmlAnnot() throws IOException {
	final File annot =super.createTmpFile(".xml");
	final PrintWriter pw = new PrintWriter(annot);
	pw.println("<x><y ncbi-gene-id=\"6331\">Hello</y></x>");
	pw.flush();
	pw.close();
	super.assertIsXml(annot);
	final File out =super.createTmpFile(".xml");
	Assert.assertEquals(new NcbiGeneDump().instanceMain(new String[] {
			"-o",out.getPath(),
			"-C",annot.getPath(),
			"SCN5A","NOTCH2"
			}),0);
	Assert.assertTrue(Files.lines(out.toPath()).anyMatch(L->L.contains("Hello")));
	super.assertIsXml(out);
	}

@Test
public void testTxtAnnot() throws IOException {
	final File annot =super.createTmpFile(".txt");
	final PrintWriter pw = new PrintWriter(annot);
	pw.println("## NCBI Gene ID. 6331 SCN5A\nHello.See also SCN10A");
	pw.flush();
	pw.close();
	super.assertIsNotEmpty(annot);
	final File out =super.createTmpFile(".xml");
	Assert.assertEquals(new NcbiGeneDump().instanceMain(new String[] {
			"-o",out.getPath(),
			"-C",annot.getPath(),
			"SCN5A","NOTCH2"
			}),0);
	Assert.assertTrue(Files.lines(out.toPath()).anyMatch(L->L.contains("Hello")));
	super.assertIsXml(out);
	}

}
