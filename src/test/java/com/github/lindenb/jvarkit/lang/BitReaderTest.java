package com.github.lindenb.jvarkit.lang;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.BitSet;
import java.util.Random;

import org.testng.Assert;
import org.testng.annotations.Test;


public class BitReaderTest {
	@Test
	public void testBitReader() throws IOException {
		final Random rnd=new Random(0L);
		final int nbits=1000;
		final BitSet bitset = new BitSet(nbits);
		for(int i=0;i< nbits;i++) {
			bitset.set(i, rnd.nextBoolean());
			}
		BitReader in=new BitReader(new ByteArrayInputStream(bitset.toByteArray()));
		for(int i=0;i< nbits;i++) {
			int b= in.read();
			Assert.assertNotEquals(b, -1);
			Assert.assertTrue(b==0 || b==1);
			Assert.assertEquals((b==1),bitset.get(i),"at("+i+")/"+bitset.size()+" expect "+bitset.get(i)+" and get "+b+"="+(b==1));
			}
		Assert.assertEquals(in.read(),-1);
		}
}
