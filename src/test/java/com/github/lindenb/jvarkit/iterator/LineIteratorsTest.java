package com.github.lindenb.jvarkit.iterator;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.iterator.LineIterators;

import htsjdk.tribble.readers.LineIterator;

public class LineIteratorsTest {
	
private void check(List<String> L) {
	Assert.assertEquals(L.size(), 3);
	Assert.assertEquals(L.get(0),"A");
	Assert.assertEquals(L.get(1),"B");
	Assert.assertEquals(L.get(2),"C");
	}
private List<String> toList(LineIterator r) {
	List<String> L = new ArrayList<>();
	while(r.hasNext()) L.add(r.next());
	return L;
	}
@Test
public void testArray() {
	LineIterator li  =  LineIterators.of(new String[] {"A","B","C"});
	check(toList(li));
	}

@Test
public void testList() {
	LineIterator li  = LineIterators.of(Arrays.asList("A","B","C"));
	check(toList(li));
	}

@Test
public void testReader() throws IOException{
	LineIterator li  = LineIterators.of(new StringReader("A\nB\nC"));
	check(toList(li));
	}
}
