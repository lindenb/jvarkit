package com.github.lindenb.jvarkit.lang;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.RuntimeIOException;

public class BitNumReaderTest {
	private BitNumReader prepareTest8Bits() throws IOException{
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		try(BinaryCodec codec = new BinaryCodec(baos)) {
			for(short i=0;i< 256;i++) {
				codec.writeUByte(i);
				}
			}
		byte[] array = baos.toByteArray();
		Assert.assertEquals(array.length,256);
		final BitReader in=new BitReader(new ByteArrayInputStream(array));
		final BitNumReader bnr = new BitNumReader(in, 8);
		Assert.assertEquals(bnr.getDefaultNBits(), 8);
		Assert.assertEquals(bnr.getDefaultMax(), 255.0);
		Assert.assertTrue(bnr.getDefaultPrecision()>0);
		return bnr;
		}
	
	private void test8BitsInt(ToIntFunction<BitNumReader> extractor) throws IOException {
		BitNumReader bnr= prepareTest8Bits();
		for(int i=0;i< 256;i++) {
			final int v=extractor.applyAsInt(bnr);
			Assert.assertTrue(v>=0);
			Assert.assertEquals(v,i,"i="+i+" v="+v);
			}
		}
	private void test8BitsFloat(ToDoubleFunction<BitNumReader> extractor) throws IOException {
		BitNumReader bnr= prepareTest8Bits();
		for(int i=0;i< 256;i++) {
			final double v=extractor.applyAsDouble(bnr);
			Assert.assertTrue(v>=0);
			Assert.assertTrue(v<=1);
			Assert.assertEquals(v,i/255.0,"i="+i+" v="+v);
			}
		}
	
@Test
public void test8BitsIntDefault() throws IOException {
	test8BitsInt(BNR->{
		try {
			return BNR.nextInt();
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		});
	}

@Test
public void test8BitsIntCustom() throws IOException {
	test8BitsInt(BNR->{
		try {
			return BNR.nextInt(8);
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		});
	}

@Test
public void test8BitsFloatDefault() throws IOException {
	test8BitsFloat(BNR->{
		try {
			return BNR.nextDouble();
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		});
	}

@Test
public void test8BitsFloatCustom() throws IOException {
	test8BitsFloat(BNR->{
		try {
			return BNR.nextDouble(8);
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		});
	}

}
