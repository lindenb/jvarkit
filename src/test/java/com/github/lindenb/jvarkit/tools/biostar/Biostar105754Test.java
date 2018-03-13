package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.broad.igv.bbfile.BBFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.CloserUtil;

public class Biostar105754Test extends TestUtils {
	@Test
	public void test01() throws IOException {
		final String bwf = SRC_TEST_RESOURCE+"/Uniqueness35bp.bigWig";
		BBFileReader bbf=new BBFileReader(bwf);
		
		final File datain = createTmpFile(".txt");
		final PrintWriter pw =new PrintWriter(datain);
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
		
		final File tsvout = createTmpFile(".txt");
		Assert.assertEquals(new Biostar105754().instanceMain(new String[] {
			"-o",tsvout.getPath(),
			"-B",bwf,
			datain.getPath()
			}),0);
		assertTsvTableIsConsitent(tsvout, null);
	}
}
