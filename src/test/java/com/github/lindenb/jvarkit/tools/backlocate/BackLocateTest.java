package com.github.lindenb.jvarkit.tools.backlocate;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Optional;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReaderTest;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest({GtfReaderTest.class,LauncherTest.class})
public class BackLocateTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException {
		try {
			Optional<Path> ref = support.getGRCh37Path();
			if(!ref.isPresent()) return ;
			Path in=support.createTmpPath(".txt");
			Path out=support.createTmpPath(".txt");
			Files.write(in,"SCN5A\tT10P\n".getBytes());
			support.assertIsNotEmpty(in);		

			Assert.assertEquals(new BackLocate().instanceMain(
					new String[]{
					"-R",ref.get().toString(),
					"-o",out.toString(),
					"--gtf",support.resource("Homo_sapiens.GRCh37.87.gtf.gz"),
					in.toString(),
					}),0);
			
			support.assertIsNotEmpty(out);
			}
		finally
			{
			support.removeTmpFiles();
			}
	}
}
