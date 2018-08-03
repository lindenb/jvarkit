package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class NaiveCnvDetectorTest extends TestUtils {
@Test
public void testO1() throws IOException {
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
}
