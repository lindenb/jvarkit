package com.github.lindenb.jvarkit.tools.onekgenomes;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.IOUtil;

public class OneThousandBamTest  extends TestUtils {
	public void test01() throws IOException
	{
	final File tmpDir  = IOUtil.createTempDir("tmp", ".TMP");
		try {
		final File output = super.createTmpFile(".bam");
		final File bed = super.createTmpFile(".bed");
		final File urlList = super.createTmpFile(".list");
		final PrintWriter pw = new PrintWriter(urlList);
		pw.println("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG03078/alignment/HG03078.mapped.ILLUMINA.bwa.MSL.low_coverage.20130415.bam");
		pw.println("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12004/alignment/NA12004.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam");
		pw.flush();
		pw.close();
		IOUtils.cat("1\t10001\t10101\n1\t20001\t20101\n", bed,false);
		assertIsBed(output);
		
		Assert.assertEquals(new OneThousandBam().instanceMain(
	     		newCmd().add(
	     		"-o",output.getPath(),
	     		"-C",tmpDir,
	     		"-T",tmpDir,
	     		"-B",bed,
	     		urlList
	     		).make()),0);
		Assert.assertIsBam(output);
		
		} 
	finally
		{
		IOUtil.deleteDirectoryTree(tmpDir);
		}
	}
}
