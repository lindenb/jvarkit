package com.github.lindenb.jvarkit.tools.structvar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class IndexCovToVcfTest {

	private final TestSupport support = new TestSupport();

	
private Path createDataFile(int nsamples) throws IOException
	{
	final Path dataFile = support.createTmpPath(".tsv");
	PrintWriter pw=new PrintWriter(Files.newBufferedWriter(dataFile));
	pw.print("#chrom\tstart\tend");
	for(int i=0;i<nsamples;i++) pw.print("\tS"+i);
	pw.println();
	for(int y=0;y< 1000;++y)
		{
		pw.print("chr1\t"+(y*1000+1)+"\t"+((y+1)*1000));
		for(int i=0;i<nsamples;i++) pw.print("\t"+support.random.nextDouble()*2.0);
		pw.println();
		}
	pw.flush();
	pw.close();
	support.assertTsvTableIsConsitent(dataFile, null);
	return dataFile;
	}
	
@Test
public void test01() throws IOException
	{
	try {
		final Path dataFile=createDataFile(10);
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(new IndexCovToVcf().instanceMain(new String[] {
			"-o",out.toString(),
			dataFile.toString()}),0);
		support.assertIsVcf(out);
		}
	catch(final Throwable err) {
		support.removeTmpFiles();
		}
	}

@Test
public void testWithPed() throws IOException
	{
	try {
		int nSamples = 10;
		final Path dataFile=createDataFile(nSamples);
		final Path ped=support.createTmpPath(".ped");
		final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(ped));
		for(int i=0;i< nSamples;++i) {
			pw.println("F\tS"+i+"\t0\t0\t0\t"+(support.random.nextBoolean()?1:2));	
			}
		pw.flush();
		pw.close();
		
		
		final Path out = support.createTmpPath(".vcf");
		Assert.assertEquals(new IndexCovToVcf().instanceMain(new String[] {
			"-o",out.toString(),
			"--pedigree",ped.toString(),
			dataFile.toString()}),0);
		support.assertIsVcf(out);
		
		}
	catch(final Throwable err) {
		support.removeTmpFiles();
		}
	}


}
