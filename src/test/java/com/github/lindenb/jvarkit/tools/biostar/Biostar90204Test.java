package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.IOUtil;

public class Biostar90204Test extends TestUtils{
	@Test(dataProvider="all-sam-or-bam-files")
	public void test01(final String bam) throws IOException {
		final File manifest = createTmpFile(".mft");
		Assert.assertEquals(
			new Biostar90204().instanceMain(newCmd().
			add("--manifest").add(manifest).
			add("-n").add(100).
			add(bam).
			make()
			),0);
		assertTsvTableIsConsitent(manifest, null);
		IOUtil.slurpLines(manifest).
			stream().
			map(L->L.split("[\t]")[0]).
			map(F->new File(F)).forEach(F->{
				try {
					assertIsValidBam(F);
					Assert.assertTrue(wc(F)<=100L);
					F.delete();
					}
				catch(final IOException err) {
					 Assert.fail(F.getPath(), err);
				}
			});
		}
	}
