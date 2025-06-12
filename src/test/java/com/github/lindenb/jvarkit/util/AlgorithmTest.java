package com.github.lindenb.jvarkit.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.testng.Assert;
import org.testng.annotations.Test;

public class AlgorithmTest {
	@Test
	public void testLowerBound() {
		final Algorithm<Integer, Integer> algo = new Algorithm<>(I->I);
		final List<Integer> L= Arrays.asList(1,1,2,5,5,5,8);
		int lb = algo.lower_bound(L, 2);
		Assert.assertEquals(lb, 2);
		lb = algo.lower_bound(L, 1);
		Assert.assertEquals(lb,0);
		lb = algo.lower_bound(L,0);
		Assert.assertEquals(lb,0);
		lb = algo.lower_bound(L,4);
		Assert.assertEquals(lb,3);
		lb =algo.lower_bound(L,8);
		Assert.assertEquals(lb,6);
		lb = algo.lower_bound(L,10);
		Assert.assertEquals(lb,L.size());
		}
	
    @Test
    public void test2() {
		final Algorithm<Integer, Integer> algo = new Algorithm<>(I->I);
    	final Random rand = new Random(0L);
		final List<Integer> data = new ArrayList<>();
		for(int i=0;i< 50;i++)
			{
			int c = 5+rand.nextInt(20);
			while(c>0) { data.add(i);--c;}
			}
		final int x0 = data.indexOf(40);
		final int x1 = data.indexOf(41);
		Assert.assertTrue(x0>0);
		Assert.assertTrue(x1>x0);
		Assert.assertEquals(x0,algo.lower_bound(data,40));
		Assert.assertEquals(x1,algo.upper_bound(data,40));
		Assert.assertEquals(x0,algo.equal_range(data,40)[0]);
		Assert.assertEquals(x1,algo.equal_range(data,40)[1]);
		Assert.assertEquals(0,algo.lower_bound(data,-1));
		Assert.assertEquals(data.size(),algo.upper_bound(data,60));
    	}
    
    
	}
