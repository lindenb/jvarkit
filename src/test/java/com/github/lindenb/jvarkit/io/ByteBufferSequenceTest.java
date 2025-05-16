package com.github.lindenb.jvarkit.io;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.function.IntUnaryOperator;

public class ByteBufferSequenceTest {

    @Test
    public void testDefaultConstructor() {
        ByteBufferSequence bbs = new ByteBufferSequence();
        Assert.assertEquals(bbs.size(), 0);
    }

    @Test
    public void testConstructorWithCapacity() {
        ByteBufferSequence bbs = new ByteBufferSequence(10);
        Assert.assertEquals(bbs.size(), 0);
    }

    @Test
    public void testConstructorWithCharSequence() {
        ByteBufferSequence bbs = new ByteBufferSequence("abc");
        Assert.assertEquals(bbs.size(), 3);
        Assert.assertEquals(bbs.toString(), "abc");
    }

    @Test
    public void testConstructorWithByteArray() {
        byte[] data = {1, 2, 3};
        ByteBufferSequence bbs = new ByteBufferSequence(data);
        Assert.assertEquals(bbs.size(), 3);
        Assert.assertEquals(bbs.toByteArray(), data);
    }

    @Test
    public void testClear() {
        ByteBufferSequence bbs = new ByteBufferSequence("hello");
        bbs.clear();
        Assert.assertEquals(bbs.size(), 0);
    }

    @Test
    public void testClone() {
        ByteBufferSequence bbs = new ByteBufferSequence("abc");
        ByteBufferSequence clone = bbs.clone();
        Assert.assertEquals(clone.size(), bbs.size());
        Assert.assertEquals(clone.toString(), bbs.toString());
        Assert.assertNotSame(clone, bbs);
    }

    @Test
    public void testWriteSingleByte() {
        ByteBufferSequence bbs = new ByteBufferSequence();
        bbs.write('a');
        Assert.assertEquals(bbs.size(), 1);
        Assert.assertEquals(bbs.toString(), "a");
    }

    @Test
    public void testWriteByteArray() throws IOException {
        ByteBufferSequence bbs = new ByteBufferSequence();
        bbs.write(new byte[]{1, 2, 3}, 0, 3);
        Assert.assertEquals(bbs.size(), 3);
        Assert.assertEquals(bbs.toByteArray(), new byte[]{1, 2, 3});
    }

    @Test
    public void testEnsureCapacity() {
        ByteBufferSequence bbs = new ByteBufferSequence(2);
        bbs.write('a');
        bbs.write('b');
        bbs.write('c'); // triggers capacity increase
        Assert.assertEquals(bbs.size(), 3);
    }

    @Test
    public void testTransform() {
        ByteBufferSequence bbs = new ByteBufferSequence("abc");
        bbs.transform(c -> c + 1);
        Assert.assertEquals(bbs.toString(), "bcd");
    }

    @Test
    public void testCharAt() {
        ByteBufferSequence bbs = new ByteBufferSequence("hello");
        Assert.assertEquals(bbs.charAt(1), 'e');
    }

    @Test
    public void testSetByteAt() {
        ByteBufferSequence bbs = new ByteBufferSequence("abc");
        bbs.setByteAt(1, (byte) 'x');
        Assert.assertEquals(bbs.toString(), "axc");
    }

    @Test
    public void testSetCharAt() {
        ByteBufferSequence bbs = new ByteBufferSequence("abc");
        bbs.setCharAt(1, 'x');
        Assert.assertEquals(bbs.toString(), "axc");
    }

    @Test
    public void testSubSequence() {
        ByteBufferSequence bbs = new ByteBufferSequence("abcdef");
        CharSequence subSeq = bbs.subSequence(1, 4);
        Assert.assertEquals(subSeq.toString(), "bcd");
    }

    @Test
    public void testAppendCharSequence() throws IOException {
        ByteBufferSequence bbs = new ByteBufferSequence("hello");
        bbs.append(" world");
        Assert.assertEquals(bbs.toString(), "hello world");
    }

    @Test
    public void testStartsWith() {
        ByteBufferSequence bbs = new ByteBufferSequence("abcdef");
        Assert.assertTrue(bbs.startsWith("abc"));
        Assert.assertFalse(bbs.startsWith("xyz"));
    }

    @Test
    public void testEquals() {
        ByteBufferSequence bbs1 = new ByteBufferSequence("abc");
        ByteBufferSequence bbs2 = new ByteBufferSequence("abc");
        ByteBufferSequence bbs3 = new ByteBufferSequence("xyz");
        Assert.assertEquals(bbs1, bbs2);
        Assert.assertNotEquals(bbs1, bbs3);
    }

    @Test
    public void testHashCode() {
        ByteBufferSequence bbs1 = new ByteBufferSequence("abc");
        ByteBufferSequence bbs2 = new ByteBufferSequence("abc");
        Assert.assertEquals(bbs1.hashCode(), bbs2.hashCode());
    }

    @Test
    public void testCompareTo() {
        ByteBufferSequence bbs1 = new ByteBufferSequence("abc");
        ByteBufferSequence bbs2 = new ByteBufferSequence("abd");
        ByteBufferSequence bbs3 = new ByteBufferSequence("abcde");

        Assert.assertTrue(bbs1.compareTo(bbs2) < 0);
        Assert.assertTrue(bbs2.compareTo(bbs1) > 0);
        Assert.assertTrue(bbs1.compareTo(bbs3) < 0);
    }

    @Test
    public void testToString() {
        ByteBufferSequence bbs = new ByteBufferSequence("hello");
        Assert.assertEquals(bbs.toString(), "hello");
    }
}
