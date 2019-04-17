package com.github.lindenb.jvarkit.tools.biostar;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class Biostar90204Test {
	
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return support.toArrayArray(support.
				allSamOrBams().
				map(F->new Object[] {F}))
				;
		}
	
	@Test(dataProvider="src1")
	public void test01(final String bam) throws IOException {
		try {
		final Path manifest = support.createTmpPath(".mft");
		Assert.assertEquals(
			new Biostar90204().instanceMain(new String[] {
			"--manifest",manifest.toString(),
			"-n","100",
			bam.toString()
			} ),0);
		support.assertTsvTableIsConsitent(manifest, null);
		Files.lines(manifest).
			map(L->L.split("[\t]")[0]).
			map(F->Paths.get(F)).
			forEach(F->{
				try {
					support.assertIsValidBam(F);
					Assert.assertTrue(support.wc(F)<=100L);
					Files.deleteIfExists(F);
					}
				catch(final IOException err) {
					 Assert.fail(F.toString(), err);
				}
			});
		} 
	finally {
		support.removeTmpFiles();
		}
	}
	}
