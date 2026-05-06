package com.github.lindenb.jvarkit.tools.vcf4igv;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

import htsjdk.samtools.util.IOUtil;

public class VcfForIGVTest {


	@Test
	public void simpleTest() 
		throws IOException
		{
		final TestSupport support = new TestSupport();
		final Path dir=null;
		try {
			dir =support.createTmpDirectory();
			final Path samplesheet = support.createTmpPath("samplesheet.csv");
			try(Writer pw = Files.newBufferedWriter(samplesheet)) {
				pw.write("sample,bam,status\n");
				pw.write("S1,"+support.resource("S1.bam")+",case");
				pw.write("S2,"+support.resource("S2.bam")+",control");
				pw.write(","+support.resource("S3.bam")+",case");
				pw.write(","+support.resource("S4.bam")+",control");
				pw.flush();
				}

			
		
			final Path out = support.createTmpPath(".json");
			final VcfForIGV cmd =new VcfForIGV();
			Assert.assertEquals(0,cmd.instanceMain(new String[] {
					"--session-dir",dir.toString(),
					"-R",support.resource("rotavirus_rf.fa"),
					"--samplesheet",samplesheet.toString(),
					"-o",out.toString(),
					support.resource("rotavirus_rf.vcf.gz")
				}));
			support.assertIsNotEmpty(out);
			}
		finally {
			support.removeTmpFiles();
			if(dir!=null) IOUtil.deleteDirectoryTree(dir.toFile());
			}
		}
	