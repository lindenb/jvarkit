package com.github.lindenb.jvarkit.util.iterator;

import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

public class FilterIteratorTest {
@Test
public void test1() {
	List<Integer> array = Arrays.asList(1,10,2,2,4,5,2);
	FilterIterator<Integer> iter = new FilterIterator<>(array.iterator(),(I)->I%2!=0);
	
	Assert.assertTrue(iter.hasNext());
	int i= iter.next();
	Assert.assertEquals(i, 1);
	Assert.assertTrue(iter.hasNext());
	i = iter.next();
	Assert.assertEquals(i, 5);
	Assert.assertFalse(iter.hasNext());
	iter.close();
	}
}
