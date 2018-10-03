package com.github.lindenb.jvarkit.tools.structvar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class MergeCnvNatorTest  extends TestUtils {
	@Test
	public void test01() throws IOException {
		final int n_samples = 5 + super.random.nextInt(5);
	final List<File> inputFiles = new ArrayList<>(n_samples);
	for(int i=0;i< n_samples;++i)
		{
		final File tmp=super.createTmpFile(".tsv");
		inputFiles.add(tmp);
		final PrintWriter pw = new PrintWriter(tmp);
		int start = 63001+random.nextInt(100);
		int end = 96000+random.nextInt(100);
		
		pw.println("deletion\tchr8:"+     start+"-"+end+"\t"+(end-start)+"\t0.65778\t4.82947e-12\t1.91606e-06\t1.02821e-11\t1.59274e-05\t0.699491");
		
		 start = 2196001+random.nextInt(100);
		 end = 2285000+random.nextInt(100);

		pw.println("duplication\tchr8:" + start+"-"+end+"\t"+(end-start)+"\t2.74121e-08\t1.80524e-62\t4.64697e-08\t1.09051e-55\t0.0191489");
		pw.flush();
		pw.close();
		super.assertTsvTableIsConsitent(tmp, null);
		}
	final File out= super.createTmpFile(".vcf");
	Assert.assertEquals(new  MergeCnvNator().instanceMain(super.newCmd().add(
			"-o",out.getPath()
			 ).addAll(inputFiles).make()),0);
	
	assertIsVcf(out);
	}
}
