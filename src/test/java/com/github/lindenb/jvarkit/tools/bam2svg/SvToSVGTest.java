package com.github.lindenb.jvarkit.tools.bam2svg;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class SvToSVGTest{
private final TestSupport support = new TestSupport();

@Test
public void test01() throws IOException{
	try {
		final Path svg = support.createTmpPath(".svg");
		
		final List<String> args= new ArrayList<>();
		args.add("-R");
		args.add(support.resource("rotavirus_rf.fa"));
		args.add("--region");
		args.add("RF01:100-200");
		args.add("-o");
		args.add(svg.toString());
		for(int i=1;i<=5;i++) args.add(support.resource("S"+i+".bam"));
		
		Assert.assertEquals(new SvToSVG().instanceMain(args),0);
		}
	finally {
		support.removeTmpFiles();
		}
	}
}
