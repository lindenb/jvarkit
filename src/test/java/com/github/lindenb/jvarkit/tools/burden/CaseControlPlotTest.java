package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class CaseControlPlotTest extends TestUtils {

@Test
public void test01() throws IOException {
	final String inputFile = SRC_TEST_RESOURCE+"/rotavirus_rf.vcf.gz";
	final File xmlFile = super.createTmpFile(".xml");
	final PrintWriter pw = new PrintWriter(xmlFile);
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
	assertIsXml(xmlFile);
	
	final File ped = super.createRandomPedigreeFromFile(inputFile);

	
	final File output = super.createTmpFile(".zip");
	Assert.assertEquals(new CaseControlPlot().instanceMain(
    		newCmd().add(
    		"-o",output,
    		"--pedigree",ped,
    		"--config",xmlFile,
    		inputFile
    		).make()
    	),0);
	
	assertZip(output);
	}
	
}
