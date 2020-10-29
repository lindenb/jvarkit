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
}
