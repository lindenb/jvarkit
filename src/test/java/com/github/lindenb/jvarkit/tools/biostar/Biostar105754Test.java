package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Random;

import org.broad.igv.bbfile.BBFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.CloserUtil;

public class Biostar105754Test {
	
	private final TestSupport support = new TestSupport();
	
	@Test
	public void test01() throws IOException {
		try {
		final Random random = new Random();
		final String bwf = support.resource("Uniqueness35bp.bigWig");
		BBFileReader bbf=new BBFileReader(bwf);
		
		final Path datain = support.createTmpPath(".txt");
		final PrintWriter pw =new PrintWriter(Files.newOutputStream(datain));
		for(final String contig:bbf.getChromosomeNames())
			{
			for(int i=0;i<100;i++)
				{
				int n=random.nextInt(100);
				int l= random.nextInt(200);
				pw.println(contig+"\t"+n+"\t"+(n+l));
				}
			}
		pw.flush();
		pw.close();
		CloserUtil.close(bbf.getBBFis());
		
		final Path tsvout = support.createTmpPath(".txt");
		Assert.assertEquals(new Biostar105754().instanceMain(new String[] {
			"-o",tsvout.toString(),
			"-B",bwf,
			datain.toString()
			}),0);
		support.assertTsvTableIsConsitent(tsvout, null);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}
}
