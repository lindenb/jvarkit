package com.github.lindenb.jvarkit.tools.gtf2xml;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.IOUtil;

public class Gtf2XmlTest {
	
	@Test
	public void runSimple() throws IOException {
		File tmp=null;
		final TestSupport support = new TestSupport();
		try {
			/*
			Assert.assertEquals(new Gtf2Xml().instanceMain(new String[] {
					"-o",tmp.toString(),
					"-m",method,
					gtf
					}),0);
*/
			}
		finally {
			support.removeTmpFiles();
			}
	}
}
