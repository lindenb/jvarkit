package com.github.lindenb.jvarkit.util.iterator;

import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

public class LineIteratorTest {
	
private void check(List<String> L) {
	Assert.assertEquals(L.size(), 3);
	Assert.assertEquals(L.get(0),"A");
	Assert.assertEquals(L.get(1),"B");
	Assert.assertEquals(L.get(2),"C");
	}
	
@Test
public void testArray() {
	LineIterator li  =new LineIterator(new String[] {"A","B","C"});
	check(li.toList());
	li.close();
	}

@Test
public void testList() {
	LineIterator li  =new LineIterator(Arrays.asList("A","B","C"));
	check(li.toList());
	li.close();
	}

@Test
public void testReader() throws IOException{
	LineIterator li  =new LineIterator(new StringReader("A\nB\nC"));
	check(li.toList());
	li.close();
	}



}
