package com.github.lindenb.jvarkit.tools.samslop;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class SamSlopTest {

	private void testSamSlop(String ...args0) throws Exception {
		final TestSupport support=new TestSupport();
		try {
			final Path samout= support.createTmpPath(".sam");
			final List<String> args = new ArrayList<>();
			args.addAll(Arrays.asList(args0));
			args.add("-o"); args.add(samout.toString());
			args.add("-R"); args.add(support.resource("rotavirus_rf.fa"));
			args.add(support.resource("S1.bam"));
			Assert.assertEquals(new SamSlop().instanceMain(args.toArray(new String[args.size()]))
					,0);
			
			support.assertIsValidBam(samout);
			}
		catch(Throwable err) {
			Assert.fail("SamSlop", err);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	
	@Test
	public void test01() throws Exception {
		testSamSlop("-m","0","-M","0");
		}
	@Test
	public void test02() throws Exception {
		testSamSlop("-m","10","-M","10");
		}
	@Test
	public void test03() throws Exception {
		testSamSlop("-m","10","-M","10","--rmClip");
		}

}
