package com.github.lindenb.jvarkit.tools.onekgenomes;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import org.testng.Assert;

import htsjdk.samtools.util.IOUtil;
import org.testng.annotations.Test;

public class OneThousandBamTest{
	private final TestSupport support = new TestSupport();

	@Test
	public void test01() throws IOException
		{
		final File tmpDir  = IOUtil.createTempDir("tmp.", ".TMP");
		try {
			final Path output = support.createTmpPath(".bam");
			final Path bed = support.createTmpPath(".bed");
			final Path urlList = support.createTmpPath(".list");
			final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(urlList));
			pw.println("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03078/alignment/HG03078.mapped.ILLUMINA.bwa.MSL.low_coverage.20130415.bam");
			pw.println("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12004/alignment/NA12004.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam");
			pw.flush();
			pw.close();
			IOUtils.cat("1\t10001\t10101\n1\t20001\t20101\n", bed,false);
			support.assertIsBed(output);
			
			Assert.assertEquals(new OneThousandBam().instanceMain(new String[] {
		     		"-o",output.toString(),
		     		"-C",tmpDir.toString(),
		     		"-T",tmpDir.toString(),
		     		"-B",bed.toString(),
		     		urlList.toString()
					}),0);
			support.assertIsValidBam(output);
			} 
		finally
			{
			IOUtil.deleteDirectoryTree(tmpDir);
			support.removeTmpFiles();
			}
		}
}
