package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class NcbiGeneDumpTest {
	
private  final TestSupport support = new TestSupport();

	
@Test
public void test01() throws IOException {
	try {
	final Path out =support.createTmpPath(".xml");
	Assert.assertEquals(new NcbiGeneDump().instanceMain(new String[] {
			"-o",out.toString(),
			"SCN5A","NOTCH2"
			}),0);
	support.assertIsXml(out);
	} finally {
		support.removeTmpFiles();
	}
	}

@Test
public void testXmlAnnot() throws IOException {
	try {
	final Path annot =support.createTmpPath(".xml");
	IOUtils.cat("<x><y ncbi-gene-id=\"6331\">Hello</y></x>", annot, false);
	support.assertIsXml(annot);
	final Path out = support.createTmpPath(".xml");
	Assert.assertEquals(new NcbiGeneDump().instanceMain(new String[] {
			"-o",out.toString(),
			"-C",annot.toString(),
			"SCN5A","NOTCH2"
			}),0);
	Assert.assertTrue(Files.lines(out).anyMatch(L->L.contains("Hello")));
	support.assertIsXml(out);
	} finally {
		support.removeTmpFiles();
	}
	}

@Test
public void testTxtAnnot() throws IOException {
	try {
	final Path annot =support.createTmpPath(".txt");
	IOUtils.cat("## NCBI Gene ID. 6331 SCN5A\nHello.See also SCN10A",annot,false);
	support.assertIsNotEmpty(annot);
	final Path out = support.createTmpPath(".xml");
	Assert.assertEquals(new NcbiGeneDump().instanceMain(new String[] {
			"-o",out.toString(),
			"-C",annot.toString(),
			"SCN5A","NOTCH2"
			}),0);
	Assert.assertTrue(Files.lines(out).anyMatch(L->L.contains("Hello")));
	support.assertIsXml(out);
	} finally {
		support.removeTmpFiles();
	}
	}

}
