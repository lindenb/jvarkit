package com.github.lindenb.jvarkit.lang;


import org.testng.Assert;
import org.testng.annotations.Test;




public class StringUtilsTest {
	
	
	
@Test
public void testCompactDNA( ) {
	Assert.assertEquals(StringUtils.compactDNA("AAAAAAT"),"6AT");
	Assert.assertEquals(StringUtils.compactDNA("AA"),"AA");
	Assert.assertEquals(StringUtils.compactDNA("A"),"A");
	Assert.assertEquals(StringUtils.compactDNA("ATTT"),"A3T");
	Assert.assertEquals(StringUtils.compactDNA("ATTTC"),"A3TC");
	}
	
@Test
public void testUncompactDNA( ) {
	Assert.assertEquals(StringUtils.uncompactDNA("6AT"),"AAAAAAT");
	Assert.assertEquals(StringUtils.uncompactDNA("10AT"),"AAAAAAAAAAT");
	Assert.assertEquals(StringUtils.uncompactDNA("1A"),"A");
	Assert.assertEquals(StringUtils.uncompactDNA("AT"),"AT");
	}
	
@Test
public void test() {
	Assert.assertTrue(StringUtils.isBlank(null));
	Assert.assertTrue(StringUtils.isBlank(" \t "));
	Assert.assertTrue(StringUtils.isDouble("-1E10"));
	Assert.assertTrue(StringUtils.isDouble("1"));
	Assert.assertFalse(StringUtils.isDouble("1a"));

	Assert.assertFalse(StringUtils.isInteger("-1E10"));
	Assert.assertTrue(StringUtils.isInteger("1"));
	Assert.assertFalse(StringUtils.isInteger("1a"));

	Assert.assertEquals(StringUtils.left("ABCD", 3),"ABC");
	Assert.assertEquals(StringUtils.left("AB", 3),"AB");

	Assert.assertEquals(StringUtils.right("ABCD", 3),"BCD");
	Assert.assertEquals(StringUtils.right("AB", 3),"AB");

	
	Assert.assertTrue(StringUtils.endsWith("x.vcf",".vcf",".a"));
	Assert.assertFalse(StringUtils.endsWith("x.bam",".vcf",".a"));
	Assert.assertFalse(StringUtils.endsWith("x.bam"));

	
	Assert.assertEquals(StringUtils.unescapeC("a\\n\\tb"),"a\n\tb");
	}



@Test 
void xmlString() {
	Assert.assertEquals(StringUtils.escapeXML("aaaaa"),"aaaaa");
	Assert.assertEquals(StringUtils.escapeXML("a<a"),"a&lt;a");
	}
}
