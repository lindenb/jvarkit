package com.github.lindenb.jvarkit.tools.bedstats;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BedStatsTest {
	
	private Path generateBed(TestSupport support) throws IOException {
		Path bed = support.createTmpPath(".bed");
		try(BufferedWriter w=Files.newBufferedWriter(bed)) {
		for(int i=0;i< 10_000;i++) {
			int B= support.random.nextInt(1000);
			int E= B + 1+ support.random.nextInt(100);
			w.append("RF0"+(1+support.random.nextInt(5)));
			w.append("\t");
			w.append(String.valueOf(B));
			w.append("\t");
			w.append(String.valueOf(E));
			w.append("\n");
		}
		w.flush();
		}
		return bed;
	}
	@Test
	public void oneBed() throws Exception {
		TestSupport support=new TestSupport();
		try {
			Path bed = generateBed(support);
			Path out = support.createTmpPath(".txt");
			
			Assert.assertEquals(
				new BedStats().instanceMain(new String[] {
						"-o",out.toString(),
						bed.toString()
						}),0);
		} finally {
			support.removeTmpFiles();
		}
	}
	@Test
	public void twoBed() throws Exception {
		TestSupport support=new TestSupport();
		try {
			Path bed1 = generateBed(support);
			Path bed2 = generateBed(support);
			Path out = support.createTmpPath(".txt");
			
			Assert.assertEquals(
				new BedStats().instanceMain(new String[] {
						"-o",out.toString(),
						bed1.toString(),
						bed2.toString()
						}),0);
		} finally {
			support.removeTmpFiles();
		}
	}
	
	@Test
	public void threeBed() throws Exception {
		TestSupport support=new TestSupport();
		try {
			Path bed1 = generateBed(support);
			Path bed2 = generateBed(support);
			Path bed3 = generateBed(support);
			Path out = support.createTmpPath(".txt");
			
			Assert.assertEquals(
				new BedStats().instanceMain(new String[] {
						"-o",out.toString(),
						bed1.toString(),
						bed2.toString(),
						bed3.toString(),
						}),0);
		} finally {
			support.removeTmpFiles();
		}
	}
}