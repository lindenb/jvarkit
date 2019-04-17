package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class CaseControlPlotTest  {
	private final TestSupport support = new TestSupport();

@Test
public void test01() throws IOException {
	try {
		final String inputFile = support.resource("rotavirus_rf.vcf.gz");
		final Path xmlFile = support.createTmpPath(".xml");
		final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(xmlFile));
		pw.println(
			"<?xml version=\"1.0\"?><config>" +
			"<filter id=\"f\">variant.getStart()%2==0</filter>"+
			"<maf id=\"m\"/>"+
			"<handler name=\"H\">"+
			"<filter ref=\"f\"/>"+
			"<case ref=\"m\"/>"+
			"<ctrl/>"+
			"</handler>"+
			"</config>"
			);
		pw.flush();
		pw.close();
		support.assertIsXml(xmlFile);
		
		final Path ped = support.createRandomPedigreeFromFile(inputFile);
	
		
		final Path output = support.createTmpPath(".zip");
		Assert.assertEquals(new CaseControlPlot().instanceMain(
	    		new String[] {
	    		"-o",output.toString(),
	    		"--pedigree",ped.toString(),
	    		"--config",xmlFile.toString(),
	    		inputFile
	    		}),0);
		
		support.assertZip(output);
		} 
	finally {
		support.removeTmpFiles();
		}
	}
	
}
