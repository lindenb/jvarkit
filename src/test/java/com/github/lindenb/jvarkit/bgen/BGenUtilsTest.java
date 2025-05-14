package com.github.lindenb.jvarkit.bgen;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.BitSet;
import java.util.Random;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;


public class BGenUtilsTest extends TestSupport{
	
	@Test
	public void testByteBuffer() throws IOException {
		final BGenUtils.ByteBuffer bf=new BGenUtils.ByteBuffer();
		final Random rnd=new Random(0L);
		final int nbytes=1007;
		for(int i=0;i< nbytes;++i) {
			bf.write((int)'A' + rnd.nextInt(26));
			}
		Assert.assertEquals(bf.size(),nbytes);
		try(InputStream is = bf.toByteArrayInputStream()) {
			byte[] array=is.readAllBytes();
			Assert.assertEquals(array.length,nbytes);
			for(byte a:array) {
				Assert.assertTrue(Character.isLetter(a));
				Assert.assertTrue(Character.isUpperCase(a));
				}
			}

		Assert.assertEquals(bf.size(),nbytes);
		}
	
	@Test
	public void testBitReader() throws IOException {
		final Random rnd=new Random(0L);
		final int nbits=1000;
		final BitSet bitset = new BitSet(nbits);
		for(int i=0;i< nbits;i++) {
			bitset.set(i, rnd.nextBoolean());
			}
		BGenUtils.BitReader in=new BGenUtils.BitReader(new ByteArrayInputStream(bitset.toByteArray()));
		for(int i=0;i< nbits;i++) {
			int b= in.read();
			Assert.assertNotEquals(b, -1);
			Assert.assertTrue(b==0 || b==1);
			Assert.assertEquals((b==1),bitset.get(i),"at("+i+")/"+bitset.size()+" expect "+bitset.get(i)+" and get "+b+"="+(b==1));
			}
		Assert.assertEquals(in.read(),-1);
		}
	@Test
	public void testBitWriter() throws IOException  {
		 final BGenUtils.ByteBuffer bf=new BGenUtils.ByteBuffer();

		final Random rnd=new Random(0L);
		final int nbits=1000;
		final BGenUtils.BitWriter w = new BGenUtils.BitWriter(bf);
		final BitSet bitset = new BitSet(nbits);
		for(int i=0;i< nbits;i++) {
				bitset.set(i, rnd.nextBoolean());
				w.write(bitset.get(i));
				}
		w.flush();
		
			
		BGenUtils.BitReader in=new BGenUtils.BitReader(bf.toByteArrayInputStream());
		for(int i=0;i< nbits;i++) {
			int b= in.read();
			Assert.assertNotEquals(b, -1);
			Assert.assertTrue(b==0 || b==1);
			Assert.assertEquals((b==1),bitset.get(i),"at("+i+")/"+bitset.size()+" expect "+bitset.get(i)+" and get "+b+"="+(b==1));
			}
		Assert.assertEquals(in.read(),-1);
		}
	}
