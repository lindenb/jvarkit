package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;

import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class IndexCovJfxTest {
	
@Test(enabled=false)
public void test01() throws IOException
	{
	/*
	final File dataFile = super.createTmpFile(".tsv");
	int nsamples=3;
	PrintWriter pw=new PrintWriter(dataFile);
	pw.print("#chrom\tstart\tend");
	for(int i=0;i<nsamples;i++) pw.print("\tS"+i);
	pw.println();
	for(int y=0;y< 1000;++y)
		{
		pw.print("chr1\t"+(y*1000+1)+"\t"+((y+1)*1000));
		for(int i=0;i<nsamples;i++) pw.print("\t"+super.random.nextDouble()*2.0);
		pw.println();
		}
	pw.flush();
	pw.close();
	super.assertTsvTableIsConsitent(dataFile, null);
	
	testJfxApplication(IndexCovJfx.class,
    		newCmd().add(
    		"--testng",
    		dataFile
    		).make()
			);
	*/
	}
}
