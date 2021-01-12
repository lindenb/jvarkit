package com.github.lindenb.jvarkit.tools.bedtools;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class BedClusterTest {
	private final TestSupport support = new TestSupport();
	
	private Path generateBed() throws IOException {
		Path bed = support.createTmpPath(".bed");
		BufferedWriter w=Files.newBufferedWriter(bed);
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
		w.close();
		return bed;
	}


	private void execute(boolean compress,boolean dict,boolean merge,boolean interval,boolean contig) throws IOException {
		try {
			Path bed = generateBed();
			Path out = support.createTmpPath(".zip");
			final List<String> args = new ArrayList<>();
			args.add("-o");
			args.add(out.toString());
			if(compress) args.add("--compress");
			if(merge) args.add("--merge");
			if(dict) { args.add("--reference");args.add(support.resource("rotavirus_rf.fa"));}
			if(interval) args.add("--interval-list");
			if(contig) args.add("--contig");
			args.add("--jobs");
			args.add("10");
			
			args.add(bed.toString());
			
			Assert.assertEquals(
					new BedCluster().instanceMain(args),0);
			support.assertZip(out);
		} finally {
			support.removeTmpFiles();
		}
	}
	
	@Test
	public void test01() throws IOException {
		execute(true,true,true,false,true);
		execute(true,false,true,false,true);
		execute(true,true,true,false,false);
		execute(false,true,false,true,true);
		execute(false,true,false,true,false);
	}
	
	@Test
	public void test02() throws IOException {
		try {
			Path bed = generateBed();
			Path out = support.createTmpPath(".zip");
			final List<String> args = new ArrayList<>();
			args.add("--size");
			args.add("1000");
			args.add("-o");
			args.add(out.toString());
			
			args.add(bed.toString());
			
			Assert.assertEquals(
					new BedCluster().instanceMain(args),0);
			support.assertZip(out);
		} finally {
			support.removeTmpFiles();
		}
	}
}
