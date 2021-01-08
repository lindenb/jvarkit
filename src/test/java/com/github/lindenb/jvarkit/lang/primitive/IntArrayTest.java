package com.github.lindenb.jvarkit.lang.primitive;

import java.util.Arrays;

import org.testng.Assert;
import org.testng.annotations.Test;

public class IntArrayTest {

	@Test
	public void test01() {
		IntArray a = new IntArray(100);
		Assert.assertTrue(a.isEmpty());
		Assert.assertTrue(a.size()==0);
		
		a = new IntArray(new int[]{0,2,4});
		Assert.assertFalse(a.isEmpty());
		Assert.assertTrue(a.size()==3);
		a.add(6);
		a.add(8);
		Assert.assertFalse(a.isEmpty());
		Assert.assertTrue(a.size()==5);
		for(int i=0;i< 5;i++) Assert.assertTrue(a.get(i)==i*2);
		
		a= new IntArray(a);
		Assert.assertFalse(a.isEmpty());
		Assert.assertTrue(a.size()==5);
		for(int i=0;i< a.size();i++) Assert.assertTrue(a.get(i)==i*2);
		
		for(int i=0;i< a.size();i++)a.set(i,i*3);
		for(int i=0;i< a.size();i++) Assert.assertTrue(a.get(i)==i*3);

		Assert.assertTrue(Arrays.equals(a.toArray(),new int[]{0,3,6,9,12}));
		a.remove(0);
		Assert.assertFalse(a.isEmpty());
		Assert.assertTrue(a.size()==4);
		Assert.assertTrue(Arrays.equals(a.toArray(),new int[]{3,6,9,12}));
	
		a.remove(a.size()-1);
		Assert.assertFalse(a.isEmpty());
		Assert.assertTrue(a.size()==3);
		Assert.assertTrue(Arrays.equals(a.toArray(),new int[]{3,6,9}));

		a.remove(1);
		Assert.assertFalse(a.isEmpty());
		Assert.assertTrue(a.size()==2);
		Assert.assertTrue(Arrays.equals(a.toArray(),new int[]{3,9}));
		
		a.insert(1,100);
		Assert.assertFalse(a.isEmpty());
		Assert.assertTrue(a.size()==3);
		Assert.assertTrue(Arrays.equals(a.toArray(),new int[]{3,100,9}));

	}
	
	
	
}
