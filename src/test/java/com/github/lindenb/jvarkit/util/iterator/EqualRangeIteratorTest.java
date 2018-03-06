package com.github.lindenb.jvarkit.util.iterator;

import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

public class EqualRangeIteratorTest {
@Test
public void test1() {
	List<Integer> array = Arrays.asList(1,2,2,2,4);
	EqualRangeIterator<Integer> iter = new EqualRangeIterator<>(array.iterator(),Integer::compareTo);
	Assert.assertTrue(iter.hasNext());
	List<Integer>  L1= iter.next();
	Assert.assertEquals(L1.size(), 1);
	Assert.assertEquals(L1.get(0).intValue(), 1);
	
	Assert.assertTrue(iter.hasNext());
	L1= iter.next();
	Assert.assertEquals(L1.size(),3);
	Assert.assertEquals(L1.get(0).intValue(), 2);
	Assert.assertEquals(L1.get(1).intValue(), 2);
	Assert.assertEquals(L1.get(2).intValue(), 2);
		
	Assert.assertTrue(iter.hasNext());
	L1= iter.next();
	Assert.assertEquals(L1.size(), 1);
	Assert.assertEquals(L1.get(0).intValue(), 4);
	
	Assert.assertFalse(iter.hasNext());
	iter.close();
	}

@Test
public void test2() {
	List<String> array = Arrays.asList("xBaa","cBccc","dBdddd","dZaa");
	EqualRangeIterator<String> iter = new EqualRangeIterator<>(array.iterator(),(S1,S2)->((int)S1.charAt(1)-(int)S2.charAt(1)));
	Assert.assertTrue(iter.hasNext());
	List<String>  L1= iter.next();
	Assert.assertEquals(L1.size(), 3);	
	Assert.assertTrue(iter.hasNext());
	L1= iter.next();
	Assert.assertEquals(L1.size(),1);
	iter.close();
	}

@Test(expectedExceptions= {IllegalStateException.class})
public void test3() {
	List<Integer> array = Arrays.asList(100,10,1);
	EqualRangeIterator<Integer> iter = new EqualRangeIterator<>(array.iterator(),Integer::compareTo);
	iter.next();
	iter.close();
	}


@Test
public void testNextTarget() {
	List<Integer> array = Arrays.asList(1,2,2,2,4,4);
	EqualRangeIterator<Integer> iter = new EqualRangeIterator<>(array.iterator(),Integer::compareTo);
	Assert.assertTrue(iter.hasNext());
	
	List<Integer>  L1= iter.next(2);
	Assert.assertEquals(L1.size(), 3);
	Assert.assertEquals(L1.get(0).intValue(), 2);
	
	L1= iter.next(1);
	Assert.assertTrue(L1.isEmpty());

	L1= iter.next(3);
	Assert.assertTrue(L1.isEmpty());
	
	L1= iter.next(4);
	Assert.assertEquals(L1.size(),2);
	
	L1= iter.next(5);
	Assert.assertTrue(L1.isEmpty());
	
	iter.close();
	}

}
