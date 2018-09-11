package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

public class Biostar336589Test extends TestUtils {
	
private File createBed() throws IOException {
	final File bed = super.createTmpFile(".bed");
	final PrintWriter pw = new PrintWriter(bed);
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
	assertIsBed(bed);
	return bed;
}
	
@Test
public void test01() throws IOException {
	final File bed = createBed();
	
	
	final File out = super.createTmpFile(".svg");
	Assert.assertEquals(new Biostar336589().instanceMain(new String[] {
			"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.dict",
			"-o",out.getPath(),
			bed.getPath(),
			}),0);
	assertIsXml(out);
	}

@Test
public void testMultiple() throws IOException {
	
	final File bed1 =  createBed();
	final File bed2 =  createBed();
	final File bed3 =  createBed();
	
	
	final File out = super.createTmpFile(".svg");
	Assert.assertEquals(new Biostar336589().instanceMain(new String[] {
			"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.dict",
			"-o",out.getPath(),
			bed1.getPath(),
			bed2.getPath(),
			bed3.getPath(),
			}),0);
	assertIsXml(out);
	}

@Test
public void testLinear() throws IOException {
	final File bed = createBed();
	final File out = super.createTmpFile(".svg");
	Assert.assertEquals(new Biostar336589().instanceMain(new String[] {
			"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.dict",
			"--width","1000",
			"-o",out.getPath(),
			bed.getPath(),
			}),0);
	assertIsXml(out);
	}

}
