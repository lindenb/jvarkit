package com.github.lindenb.jvarkit.tools.sashimi;

import java.io.IOException;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class PlotSashimiTest {
	final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException {
		
		try {
			Path out = support.createTmpPath(".zip");
			Path mf = support.createTmpPath(".mf");
			
			Assert.assertEquals(new PlotSashimi().instanceMain(new String[] {
					"-o",out.toString(),
					"-r","chr3:38597400-38599300",
					"-m",mf.toString(),
					support.resource("ENCFF331CGL.rnaseq.b38.bam")
				}),0
				);	
			support.assertTsvTableIsConsitent(mf,null);
			support.assertZip(out);
		} finally {
			support.removeTmpFiles();
		}
	}

}
