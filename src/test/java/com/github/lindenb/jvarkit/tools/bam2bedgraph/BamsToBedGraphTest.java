package com.github.lindenb.jvarkit.tools.bam2bedgraph;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BamsToBedGraphTest {
	
	@Test
	public void testInterval() throws IOException {
		final TestSupport support = new TestSupport();
		try {
			Path bedOut = support.createTmpPath(".bed");
			
			final List<String> args= new ArrayList<>();
			args.add("--interval");
			args.add("RF01:1-1000");
			args.add("-R");
			args.add(support.resource("rotavirus_rf.fa"));
			args.add("-o");
			args.add(bedOut.toString());
			Arrays.asList("1","2","3","4","5").stream().
			map(S->support.resource("/S"+S+".bam")).
			forEach(S->args.add(S));
			
			Assert.assertEquals(new BamsToBedGraph().instanceMain(
					args.toArray(new String[args.size()])),
					0);
			support.assertIsBed(bedOut);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}


	@Test
	public void testDictOnly() throws IOException {
		final TestSupport support = new TestSupport();
		try {
			Path bedOut = support.createTmpPath(".bed");
			
			final List<String> args= new ArrayList<>();
			args.add("-R");
			args.add(support.resource("rotavirus_rf.fa"));
			args.add("-o");
			args.add(bedOut.toString());
			Arrays.asList("1","2","3","4","5").stream().
			map(S->support.resource("/S"+S+".bam")).
			forEach(S->args.add(S));
			
			Assert.assertEquals(new BamsToBedGraph().instanceMain(
					args.toArray(new String[args.size()])),
					0);
			support.assertIsBed(bedOut);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}

}