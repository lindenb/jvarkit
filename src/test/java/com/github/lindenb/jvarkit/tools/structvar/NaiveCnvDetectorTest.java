package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class NaiveCnvDetectorTest {
	
	private final TestSupport support = new TestSupport();
	
	
@Test
public void testOneFile() throws IOException {
	try {
	final Path out= support.createTmpPath(".tsv");
	final Path tmp= support.createTmpPath(".tmp");
	final PrintWriter pw=IOUtils.openPathForPrintWriter(tmp);
	for(int i=1;i< 100_000;++i)
		{
		pw.print("chr1\t");
		pw.print(i);
		for(int j=0;j< 30;++j)
			{
			int depth=50+(support.random.nextInt(10)*(support.random.nextBoolean()?1:-1));
			if(j==40 && i> 10_000 && i<20_000) depth/=2;
			pw.print("\t");
			pw.print(depth);
			}
		pw.println();
		}
	pw.flush();
	pw.close();
	Assert.assertEquals(new NaiveCnvDetector().instanceMain(new String[] {
			"-o",out.toString(),
			tmp.toString()}),
			0);
	support.assertTsvTableIsConsitent(out, null);
	} finally 
	{
		support.removeTmpFiles();
	}
	}

@Test
public void testMultipleFiles() throws IOException {
	try {
		final Path listF= support.createTmpPath(".list");
		final List<Path> tmps= new ArrayList<>(10);
		for(int x=0;x<10;++x) {
			final Path tmp = support.createTmpPath(".tmp");
			tmps.add(tmp);
			final PrintWriter pw=IOUtils.openPathForPrintWriter(tmp);
			for(int i=1;i< 100_000;++i)
				{
				pw.print("chr1\t");
				pw.print(i);
				pw.print("\t");
				int depth=50+(support.random.nextInt(10)*(support.random.nextBoolean()?1:-1));
				if(x==5 && i> 10_000 && i<20_000) depth/=2;
				pw.print(depth);
				pw.println();
				}
			pw.flush();
			pw.close();
			}
		PrintWriter pw=IOUtils.openPathForPrintWriter(listF);
		for(final Path tmpF : tmps)
			pw.println(tmpF.toString());
		pw.flush();
		pw.close();
		
		final Path dict=support.createTmpPath(".dict");
		pw=IOUtils.openPathForPrintWriter(dict);
		pw.println("@HD\tVN:1.5\tSO:unsorted");
		pw.println("@SQ\tSN:chr1\tLN:10000000");
		pw.flush();
		pw.close();
		
		
		final Path out= support.createTmpPath(".tsv");
		Assert.assertEquals(new NaiveCnvDetector().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",dict.toString(),
				listF.toString()}),
				0);
		support.assertTsvTableIsConsitent(out, null);
		}
	finally 
		{
		support.removeTmpFiles();
		}
	}


}
