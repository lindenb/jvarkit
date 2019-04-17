package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Random;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtilsTest;


@AlsoTest({VCFUtilsTest.class,BedLineCodec.class})
public class Biostar336589Test  {

private final TestSupport support = new TestSupport();
private final Random random = new Random();
	
private Path createBed() throws IOException {
	final Path bed = support.createTmpPath(".bed");
	final PrintWriter pw = new PrintWriter(Files.newOutputStream(bed));
	for(int i=0;i< 10;i++)
		{
		int start= random.nextInt(600);
		int end = start + random.nextInt(60);
		char strand;
		switch(random.nextInt(3))
			{
			case 0: strand='+'; break;
			case 1: strand='-'; break;
			default:strand='.'; break;
			}
		pw.println("RF"+ String.format("%02d", 1+random.nextInt(10))+
				"\t"+start+"\t"+end+"\trota"+random.nextInt(10000)+
				"\t"+random.nextInt(1000)+"\t"+strand+"\t.\t.\t0,255,0");
		}
	
	pw.close();
	support.assertIsBed(bed);
	return bed;
}
	
@Test
public void test01() throws IOException {
	try {
		final Path bed = createBed();
		final Path out = support.createTmpPath(".svg");
		Assert.assertEquals(new Biostar336589().instanceMain(new String[] {
				"-R",support.resource("rotavirus_rf.dict"),
				"-o",out.toString(),
				bed.toString(),
				}),0);
		support.assertIsXml(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}

@Test
public void testMultiple() throws IOException {
	try {
		Path  bed1 =  createBed();
		Path  bed2 =  createBed();
		Path  bed3 =  createBed();
		
		
		final Path out = support.createTmpPath(".svg");
		Assert.assertEquals(new Biostar336589().instanceMain(new String[] {
				"-R",support.resource("rotavirus_rf.dict"),
				"-o",out.toString(),
				bed1.toString(),
				bed2.toString(),
				bed3.toString(),
				}),0);
		support.assertIsXml(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}

@Test
public void testLinear() throws IOException {
	try {
		final Path bed = createBed();
		final Path out = support.createTmpPath(".svg");
		Assert.assertEquals(new Biostar336589().instanceMain(new String[] {
				"-R",support.resource("rotavirus_rf.dict"),
				"--width","1000",
				"-o",out.toString(),
				bed.toString()
				}),0);
		support.assertIsXml(out);
		}
	finally {
		support.removeTmpFiles();
		}
	}

}
