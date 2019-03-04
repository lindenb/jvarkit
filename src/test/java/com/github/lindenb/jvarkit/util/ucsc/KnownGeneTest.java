package com.github.lindenb.jvarkit.util.ucsc;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.CharSplitterTest;
import com.github.lindenb.jvarkit.tools.tests.AlsoTest;

@AlsoTest(CharSplitterTest.class)
public class KnownGeneTest {
	
	@Test
	public void test01() {
		String tokens[]=CharSplitter.SPACE.split("592 ENST00000379319.1 chr1 - 1018079 1041507 1018272 1026923 10 1018079,1019732,1019860,1021257,1022518,1022881,1025732,1026851,1027370,1041335, 1018367,1019763,1019886,1021392,1022584,1022977,1025808,1026945,1027483,1041507,");	
		KnownGene kg = new KnownGene(tokens);
		Assert.assertEquals(kg.getName(), "ENST00000379319.1");
		Assert.assertTrue(kg.isNegativeStrand());
		Assert.assertFalse(kg.isPositiveStrand());
		Assert.assertEquals(kg.getTxStart(),1018079);
		Assert.assertEquals(kg.getTxEnd(),1041507);
		Assert.assertEquals(kg.getCdsStart(),1018272);
		Assert.assertEquals(kg.getCdsEnd(),1026923);
		Assert.assertEquals(kg.getExonCount(),10);
		Assert.assertEquals(kg.getExonStart(0),1018079);
		Assert.assertEquals(kg.getExonEnd(0), 1018367);
		Assert.assertEquals(kg.getExon(0).getStart(),1018079);
		Assert.assertEquals(kg.getExon(0).getEnd(),1018367);
		Assert.assertEquals(kg.getIntron(0).getStart(),1018367);
		Assert.assertFalse(kg.isNonCoding());
		}
	
}
