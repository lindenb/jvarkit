package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;


import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;


public class MergeCnvNatorTest {
	
	private final TestSupport support = new TestSupport();

	
	@Test
	public void test01() throws IOException {
		try {
			final int n_samples = 5 + support.random.nextInt(5);
			final List<Path> inputFiles = new ArrayList<>(n_samples);
			for(int i=0;i< n_samples;++i)
				{
				final Path tmp=support.createTmpPath(".tsv");
				inputFiles.add(tmp);
				final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(tmp));
				int start = 63001+support.random.nextInt(100)+i;
				int end = 96000+support.random.nextInt(100)+i;
				
				pw.println("deletion\tchr8:"+     start+"-"+end+"\t"+(end-start)+"\t0.65778\t4.82947e-12\t1.91606e-06\t1.02821e-11\t1.59274e-05\t0.699491");
				
				 start = 2196001+support.random.nextInt(100);
				 end = 2285000+support.random.nextInt(100);
		
				pw.println("duplication\tchr8:" + start+"-"+end+"\t"+(end-start)+"\t0.65778\t4.82947e-12\t1.91606e-06\t1.02821e-11\t1.59274e-05\t0.699491");
				pw.flush();
				pw.close();
				support.assertTsvTableIsConsitent(tmp, null);
				}
			final Path out= support.createTmpPath(".vcf");
			
			final List<String> args= new ArrayList<String>();
			args.add("-o");
			args.add(out.toString());
			inputFiles.stream().map(P->P.toString()).forEach(S->args.add(S));
			
			Assert.assertEquals(new  MergeCnvNator().instanceMain(args),0);
			
			support.assertIsVcf(out);
			}
		finally {
			support.removeTmpFiles();
		}
	}
}
