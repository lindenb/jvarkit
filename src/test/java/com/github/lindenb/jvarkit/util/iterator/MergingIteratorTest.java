package com.github.lindenb.jvarkit.util.iterator;

import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.Test;

public class MergingIteratorTest {
@Test
public void test1() {
	final MergingIterator<Integer> iter = new MergingIterator<>(
			Integer::compare,
			Arrays.asList(
					Arrays.asList(1,3,5,7,9).iterator(),
					Arrays.asList(0,2,4,6,8,10).iterator(),
					Arrays.asList(1,3,5,7,9).iterator(),
					Arrays.asList(0,2,4,6,8,10).iterator(),
					Arrays.asList(11,11,12,12).iterator()
			));
	for(int i=0;i<=12;i++)
		{
		for(int n=0;n<2;++n) {
			Assert.assertTrue(iter.hasNext());
			int j= iter.next();
			Assert.assertEquals(i,j);
			}
		}
	Assert.assertFalse(iter.hasNext());
	iter.close();
	}
@Test(expectedExceptions= {IllegalStateException.class})
public void test2() {
	final MergingIterator<Integer> iter = new MergingIterator<>(
			Integer::compare,
			Arrays.asList(
					Arrays.asList(10,9).iterator(),
					Arrays.asList(4,2).iterator()
			));
	while(iter.hasNext()) iter.next();
	iter.close();
	}
}
