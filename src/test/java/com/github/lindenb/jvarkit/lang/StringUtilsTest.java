package com.github.lindenb.jvarkit.lang;


import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;

@AlsoTest({AbstractCharSequenceTest.class,DelegateCharSequenceTest.class})
public class StringUtilsTest {
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
}
