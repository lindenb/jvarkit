package com.github.lindenb.jvarkit.tools.viewmate;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

@AlsoTest(LauncherTest.class)
public class SamViewWithMateTest {
	private final TestSupport support =new TestSupport();


private void test01(boolean withIndex) throws IOException {
	final Path input = Paths.get(support.resource("HG02260.transloc.chr9.14.bam"));
	final Path output = support.createTmpPath(".bam");
	
	List<String> args = new ArrayList<>();
	args.add("-o");args.add(output.toString());
	args.add("-r");args.add("9:137230721-137230796");
	if(!withIndex) args.add("--streaming");
	args.add(input.toString());
	
	Assert.assertEquals(new SamViewWithMate().instanceMain(args),0);
    support.assertIsValidBam(output);
    SamReaderFactory srf=SamReaderFactory.makeDefault();
	for(final String str: new String[]{"9","14"}) {
		try(SamReader sr = srf.open(output)) {
			Assert.assertTrue(sr.iterator().stream().anyMatch(R->str.equals(R.getContig())));
			}
		}
    
	}

@Test
public void testWithIndex() throws IOException {
	try {
		test01(true);
		}
	finally
		{
			support.removeTmpFiles();
		}
	}
@Test
public void testWithoutIndex() throws IOException {
	try {
		test01(false);
		}
	finally
		{
			support.removeTmpFiles();
		}
	}

}
