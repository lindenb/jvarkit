package com.github.lindenb.jvarkit.lang;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.AlsoTest;

@AlsoTest(AbstractCharSequenceTest.class)
public class DelegateCharSequenceTest {
	@Test
	public void testCharSeq() {
		final String s1 = "0123456789";
		final DelegateCharSequence abs = new DelegateCharSequence(s1);
		Assert.assertEquals(s1.length(), abs.length());
		Assert.assertEquals(s1, abs.toString());
		for(int i=0;i< s1.length();i++) {
			Assert.assertEquals(s1.charAt(i), abs.charAt(i));	
			}
		Assert.assertTrue(abs.hasStringValue(s1));
		}
}
