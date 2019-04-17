package com.github.lindenb.jvarkit.util;

import org.testng.Assert;
import org.testng.annotations.Test;

public class CounterTest {

@Test
public void testWithInt() {
	final Counter<Integer> counter=new Counter<>();
	for(int x : new int[]{1,1,1,2,2,10,10})
		{
		counter.incr(x);
		}
	Assert.assertEquals(counter.count(1), 3);
	Assert.assertEquals(counter.count(2), 2);
	Assert.assertEquals(counter.count(10), 2);
	Assert.assertEquals(counter.count(3), 0);
	Assert.assertEquals(counter.getTotal(), 7);
	Assert.assertEquals(counter.getMostFrequent(),Integer.valueOf(1));
	}
}
