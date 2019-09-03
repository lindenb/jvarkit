package com.github.lindenb.jvarkit.samtools;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMReadGroupRecord;

public class SAMReadGroupParserTest {
@Test
public void test01() {
	final SAMReadGroupParser p=new SAMReadGroupParser();
	SAMReadGroupRecord rg= p.apply("@RG\tID:foo\tSM:bar");
	Assert.assertNotNull(rg);
	Assert.assertEquals(rg.getId(),"foo");
	Assert.assertEquals(rg.getAttribute("SM"),"bar");
	Assert.assertEquals(rg.getSample(),"bar");
	}
}
