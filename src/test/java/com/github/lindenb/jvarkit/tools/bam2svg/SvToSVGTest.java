package com.github.lindenb.jvarkit.tools.bam2svg;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import org.testng.Assert;
import org.testng.annotations.Test;
import com.github.lindenb.jvarkit.tools.tests.TestUtils;


public class SvToSVGTest extends TestUtils{
@Test
public void test01() throws IOException{
	final File svg = super.createTmpFile(".svg");
	
	Assert.assertEquals(new SvToSVG().instanceMain(newCmd().add(
			"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.fa",
			"--region","RF01:100-200",
			"-o",svg).
			add(Arrays.asList("1","2","3","4","5").stream().
				map(S->SRC_TEST_RESOURCE+"/S"+S+".bam").
				toArray()).
			make()
			),0);
	assertIsXml(svg);
	}
}
