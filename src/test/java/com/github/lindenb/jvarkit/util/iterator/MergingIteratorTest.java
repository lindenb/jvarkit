package com.github.lindenb.jvarkit.util.iterator;

import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.Test;

public class MergingIteratorTest {
@Test
public void test1() {
	MergingIterator<Integer> iter = new MergingIterator<>(
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
}
