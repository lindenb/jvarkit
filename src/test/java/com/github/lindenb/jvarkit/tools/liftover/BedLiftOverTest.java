package com.github.lindenb.jvarkit.tools.liftover;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class BedLiftOverTest {

	private Path  mockChain(Path out) throws IOException {
		Files.writeString(out, "chain 1000 RF01 3302 + 0 2600 RF02 2687 + 0 2600 1\n2600\t0\t0\n0");
		return out;
		}

	@Test
	public void simpleLift() {
		final TestSupport support =new TestSupport();
		try {
			Path chain = mockChain(support.createTmpPath(".chain"));
			Path inbed = support.createTmpPath(".bed");
			Path outbed = support.createTmpPath(".bed");
			Files.writeString(inbed, "RF01 1 100 Hello\n".replace(" ", "\t"));

			List<String> args = Arrays.asList(
					"--chain",chain.toAbsolutePath().toString(),
					"-R",support.resource("rotavirus_rf.fa"),
					"-o",outbed.toAbsolutePath().toString(),
					inbed.toAbsolutePath().toString()
					);
			
			Assert.assertEquals(new BedLiftOver().instanceMain(args),0);
			Assert.assertEquals( IOUtils.slurpPath(outbed),"RF02\t1\t100\tHello\n");
			}
		catch(Throwable err) {
			Assert.fail("boum",err);
			}
		finally {
			support.removeTmpFiles();
		}
	}
	
	@Test
	public void shiftColumns() throws IOException {
		final TestSupport support =new TestSupport();
		try {
			Path chain = mockChain(support.createTmpPath(".chain"));
			Path inbed = support.createTmpPath(".bed");
			Path outbed = support.createTmpPath(".bed");
			Files.writeString(inbed, "A RF01 1 100 Hello\n".replace(" ", "\t"));
			Assert.assertEquals(new BedLiftOver().instanceMain(Arrays.asList(
				"--chain",chain.toAbsolutePath().toString(),
				"-R",support.resource("rotavirus_rf.fa"),
				"-o",outbed.toAbsolutePath().toString(),
				"-c","2,3,4",
				inbed.toAbsolutePath().toString()
				)),0);
			Assert.assertEquals( IOUtils.slurpPath(outbed),"A\tRF02\t1\t100\tHello\n");
			}
		finally {
			support.removeTmpFiles();
		}
	}	
	
}
