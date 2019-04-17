package com.github.lindenb.jvarkit.tools.bam2svg;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class WesCnvSvgTest{
	private final TestSupport support = new TestSupport();

@Test
public void testWithBed() throws IOException{
	try {
	
	
	final Path svg = support.createTmpPath(".svg");
	final Path bed = support.createTmpPath(".bed");
	final PrintWriter pw=new PrintWriter(Files.newOutputStream(bed));
	pw.println("RF01\t1\t100");
	pw.println("RF02\t1\t100");
	pw.println("RF03\t1\t100");
	pw.flush();
	pw.close();
	support.assertIsBed(bed);
	
	final List<String> args= new ArrayList<>();
	args.add("-R");
	args.add(support.resource("rotavirus_rf.fa"));
	args.add("-B");
	args.add(bed.toString());
	args.add("-o");
	args.add(svg.toString());
	for(int i=1;i<=5;i++) args.add(support.resource("S"+i+".bam"));
	
	Assert.assertEquals(new WesCnvSvg().instanceMain(args),0);
	support.assertIsXml(svg);
	} finally {
		support.removeTmpFiles();
	}
	}
@Test
public void testWithRegion() throws IOException{
	try {
	final Path svg = support.createTmpPath(".svg");
	
	final List<String> args= new ArrayList<>();
	args.add("-R");
	args.add(support.resource("rotavirus_rf.fa"));
	args.add("--region");
	args.add("RF01:1-100;RF02:200-300");
	args.add("-o");
	args.add(svg.toString());
	for(int i=1;i<=5;i++) args.add(support.resource("S"+i+".bam"));
	
	Assert.assertEquals(new WesCnvSvg().instanceMain(args),0);
	support.assertIsXml(svg);
	} finally {
		support.removeTmpFiles();
	}
}

}
