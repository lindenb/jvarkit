package com.github.lindenb.jvarkit.math;

import org.testng.Assert;
import org.testng.annotations.Test;

public class DiscreteMedianTest {

@Test
public void test1() {
	DiscreteMedian<Integer> d = new DiscreteMedian<>();
	for(int i=0;i< 100;i++) d.add(100);
	for(int i=0;i< 10;i++) d.add(1);
	int m= (int)d.getMedian().orElse(0.0);
	Assert.assertEquals(m,100);
	}
@Test
public void test2() {
	DiscreteMedian<Integer> d = new DiscreteMedian<>();
	for(int i=0;i< 100;i++) d.add(10);
	for(int i=0;i< 100;i++) d.add(20);
	int m= (int)d.getMedian().orElse(0.0);
	Assert.assertEquals(m,15);
	}
@Test
public void test3() {
	DiscreteMedian<Integer> d = new DiscreteMedian<>();
	d.add(0);
	d.add(1);
	double m= d.getMedian().orElse(0.0);
	Assert.assertTrue(0.49 < m && m  < 0.51);
	}

@Test
public void test4() {
	DiscreteMedian<Integer> d = new DiscreteMedian<>();
	d.add(10000);
	for(int i=0;i< 5 ;++i) d.add(10);
	int m= (int)d.getMin().getAsDouble();
	Assert.assertEquals(m, 10);
	Assert.assertEquals((int)d.getMedian().orElse(0.0), 10);
	}

@Test
public void test5() {
	DiscreteMedian<Integer> d = new DiscreteMedian<>();
	for(int i=0;i< 10000;i++) d.add(10);
	d.add(199);
	int m= (int)d.getMax().getAsDouble();
	Assert.assertEquals(m, 199);
	Assert.assertEquals((int)d.getMedian().orElse(0.0), 10);
	}

public void test6() {
	DiscreteMedian<Integer> d = new DiscreteMedian<>();
	for(int i=0;i< 10000;i++) d.add(1000);
	for(int i=0;i< 10000;i++) d.add(10);
	for(int i=0;i< 10000;i++) d.add(100000);
	double m= (int)d.getMedian().orElse(0.0);
	Assert.assertEquals(m,1000);
	}

}
