package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class NaiveCnvDetectorTest extends TestUtils {
@Test
public void testOneFile() throws IOException {
	final File out= super.createTmpFile(".tsv");
	final File tmp= super.createTmpFile(".tmp");
	final PrintWriter pw=new PrintWriter(tmp);
	for(int i=1;i< 100_000;++i)
		{
		pw.print("chr1\t");
		pw.print(i);
		for(int j=0;j< 30;++j)
			{
			int depth=50+(random.nextInt(10)*(random.nextBoolean()?1:-1));
			if(j==40 && i> 10_000 && i<20_000) depth/=2;
			pw.print("\t");
			pw.print(depth);
			}
		pw.println();
		}
	pw.flush();
	pw.close();
	Assert.assertEquals(new NaiveCnvDetector().instanceMain(new String[] {
			"-o",out.getPath(),
			tmp.getPath()}),
			0);
	super.assertTsvTableIsConsitent(out, null);
	}

@Test
public void testMultipleFiles() throws IOException {
	final File listF= super.createTmpFile(".list");
	final List<File> tmps= new ArrayList<>(10);
	for(int x=0;x<10;++x) {
		final File tmp = super.createTmpFile(".tmp");
		tmps.add(tmp);
		final PrintWriter pw=new PrintWriter(tmp);
		for(int i=1;i< 100_000;++i)
			{
			pw.print("chr1\t");
			pw.print(i);
			pw.print("\t");
			int depth=50+(random.nextInt(10)*(random.nextBoolean()?1:-1));
			if(x==5 && i> 10_000 && i<20_000) depth/=2;
			pw.print(depth);
			pw.println();
			}
		pw.flush();
		pw.close();
		}
	PrintWriter pw=new PrintWriter(listF);
	for(final File tmpF : tmps)
		pw.println(tmpF.getPath());
	pw.flush();
	pw.close();
	
	final File dict=super.createTmpFile(".dict");
	pw=new PrintWriter(dict);
	pw.println("@HD\tVN:1.5\tSO:unsorted");
	pw.println("@SQ\tSN:chr1\tLN:10000000");
	pw.flush();
	pw.close();
	
	
	final File out= super.createTmpFile(".tsv");
	Assert.assertEquals(new NaiveCnvDetector().instanceMain(new String[] {
			"-o",out.getPath(),
			"-R",dict.getPath(),
			listF.getPath()}),
			0);
	super.assertTsvTableIsConsitent(out, null);
	}


}
