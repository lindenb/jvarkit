package com.github.lindenb.jvarkit.tools.sam4weblogo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

import java.nio.file.Path;
import java.nio.file.Paths;

import htsjdk.samtools.util.Interval;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;

@AlsoTest(LauncherTest.class)
public class SAM4WebLogoTest
	{
	private final TestSupport support = new TestSupport();

	
	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(
				support.allSamOrBams().
				filter(support.bamHasIndex).
				map(F->new Object[] {F})
				);
		}
	
	private Path basetest(final String inBam,String params,boolean fastq) throws IOException {
		if(SAMSequenceDictionaryExtractor.extractDictionary(Paths.get(inBam)).
				getSequences().stream().
				mapToInt(L->L.getSequenceLength()).
				max().orElse(0) > 1_000_000) return null;

		final Path out = support.createTmpPath(fastq?".fastq":".txt");
		for(final Interval interval: support.randomIntervalsFromDict(Paths.get(inBam), 20,1000)) {
			if(interval.getContig().contains(":")) continue;
			List<String> args= new ArrayList<>();
			args.add("-o");
			args.add(out.toString());
			args.add("-r");
			args.add(interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd());
			Arrays.stream(params.split("[ \t]+")).filter(S->!StringUtils.isBlank(S)).forEach(S->args.add(S));
			args.add(inBam);
			
			final int ret=new SAM4WebLogo().instanceMain(args);
			if(ret!=0) {
				Assert.fail(Arrays.asList(args).toString());
				}
			
			}
		return out;
		}
	
	@Test(dataProvider="src1")
	public void testNoClip(final String inBam) throws IOException {
		try {
			basetest(inBam,"",false);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	@Test(dataProvider="src1")
	public void testClip(final String inBam) throws IOException {
		try {
			basetest(inBam,"-c",false);
			}
		finally
			{
			support.removeTmpFiles();
			}
		}
	
	@Test(dataProvider="src1")
	public void testFastq(final String inBam) throws IOException {
		try {
			final Path out=basetest(inBam,"-c --fastq -fqu _ -fqp _",true);
			if(out!=null) support.assertIsFastq(out);
			}
			finally
				{
				support.removeTmpFiles();
				}
		}
	}
