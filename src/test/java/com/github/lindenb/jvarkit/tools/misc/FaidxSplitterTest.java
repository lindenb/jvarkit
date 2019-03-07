package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.IOUtilsTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtilsTest;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodecTest;


@AlsoTest({BedLineCodecTest.class,IOUtilsTest.class,SequenceDictionaryUtilsTest.class})
public class FaidxSplitterTest {
	
private final TestSupport support = new TestSupport();
@Test
public void test01() throws IOException {
	try {
		final Path out = support.createTmpPath(".bed");

		final Path gap = support.createTmpPath(".bed");
		IOUtils.cat("RF01\t100\t200\n",gap,false);
		support.assertIsBed(gap);
		final Path genes = support.createTmpPath(".bed");
		IOUtils.cat("G1\t10\t90\nG2\t300\t400\n",genes,false);
		support.assertIsBed(genes);
		Assert.assertEquals(new FaidxSplitter().instanceMain(new String[] {
			"-o",out.toString(),
			"-gap",gap.toString(),
			"-gene",genes.toString(),
			"-R",support.resource("rotavirus_rf.fa")
			}),0);
		support.assertIsBed(out);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}
}
