package com.github.lindenb.jvarkit.lang;

import org.testng.Assert;
import org.testng.annotations.Test;

public class AbstractCharSequenceTest {
@Test
public void testCharSeq() {
	final String s1 = "0123456789";
	final AbstractCharSequence abs = new AbstractCharSequence() {
			@Override
			public int length() {
				return s1.length();
			}
			
			@Override
			public char charAt(int index) {
				return s1.charAt(index);
			}
		};
	Assert.assertEquals(s1.length(), abs.length());
	Assert.assertEquals(s1, abs.toString());
	for(int i=0;i< s1.length();i++) {
		Assert.assertEquals(s1.charAt(i), abs.charAt(i));	
		}
	Assert.assertTrue(abs.hasStringValue(s1));
	}
}
