package com.github.lindenb.jvarkit.tools.bed2hilbert;

import java.io.IOException;
import java.io.Writer;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BedToHilbertTest {

	private Path createBed(final TestSupport support) throws IOException  {
		final Path bed = support.createTmpPath(".bed");
		try(Writer w= Files.newBufferedWriter(bed)) {
			w.write("RF01\t10\t1100\n");
			w.write("RF01\t100\t1000\n");
			w.write("RF02\t100\t500\n");
			}
		return bed;
		}

	
	@Test
	public void testWithRef()
		throws IOException
		{
		final TestSupport support = new TestSupport();
		try {
			final Path in = createBed(support);
			final Path out = support.createTmpPath(".xml");
			final BedToHilbert cmd =new BedToHilbert();
			Assert.assertEquals(0,cmd.instanceMain(new String[] {
				"-R",support.resource("rotavirus_rf.dict"),
				"-o",out.toString(),
				in.toString()
				}));
			support.assertIsXml(out);
			}
		finally {
			support.removeTmpFiles();
			}
		}
}