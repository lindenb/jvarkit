package com.github.lindenb.jvarkit.tools.bam2svg;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestUtils;

import htsjdk.samtools.util.Interval;

public class WesCnvSvgTest extends TestUtils{
@Test
public void test01() throws IOException{
	final File svg = super.createTmpFile(".svg");
	final File bed = super.createTmpFile(".bed");
	final List<Interval> interval =super.randomIntervalsFromDict(new File(SRC_TEST_RESOURCE+"/rotavirus_rf.dict"), 10);
	final PrintWriter pw=new PrintWriter(bed);
	for(Interval i:interval) {
		pw.println(i.getContig()+"\t"+(i.getStart()-1)+"\t"+i.getEnd());
	}
	pw.flush();
	pw.close();
	assertIsBed(bed);
	
	Assert.assertEquals(new WesCnvSvg().instanceMain(newCmd().add(
			"-R",SRC_TEST_RESOURCE+"/rotavirus_rf.fa",
			"-B",bed,
			"-o",svg).
			add(Arrays.asList("1","2","3","4","5").stream().
				map(S->SRC_TEST_RESOURCE+"/S"+S+".bam").
				toArray()).
			make()
			),0);
	assertIsXml(svg);
	}
}
