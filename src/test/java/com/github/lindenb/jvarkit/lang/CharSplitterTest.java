package com.github.lindenb.jvarkit.lang;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;


public class CharSplitterTest {
@Test
public void testCount() {
	
	Assert.assertEquals(CharSplitter.COMMA.countTokens(",,,,"),1);
	Assert.assertEquals(CharSplitter.COMMA.countTokens("1,2,3,4,,,,,"),4);
	Assert.assertEquals(CharSplitter.COMMA.countTokens(",,1,,2,,,,"),5);
	}
@Test
public void testSequenceList01() {
	final List<CharSequence> L=CharSplitter.COMMA.splitAsCharSequenceList(",AA,BBB,,DDDD,,");
	Assert.assertEquals(L.size(),5);
	Assert.assertEquals(L.get(0).toString(),"");
	Assert.assertEquals(L.get(1).toString(),"AA");
	Assert.assertEquals(L.get(2).toString(),"BBB");
	Assert.assertEquals(L.get(3).toString(),"");
	Assert.assertEquals(L.get(4).toString(),"DDDD");
	}
@Test
public void testSequenceList02() {
	final List<CharSequence> L=CharSplitter.COMMA.splitAsCharSequenceList(",AA,BBB,,DDDD,,",3);
	Assert.assertEquals(L.size(),3);
	Assert.assertEquals(L.get(0).toString(),"");
	Assert.assertEquals(L.get(1).toString(),"AA");
	Assert.assertEquals(L.get(2).toString(),"BBB,,DDDD");
	}
@Test
public void testSequenceList03() {
	final List<CharSequence> L=CharSplitter.COMMA.splitAsCharSequenceList("A,B",1000);
	Assert.assertEquals(L.size(),2);
	Assert.assertEquals(L.get(0).toString(),"A");
	Assert.assertEquals(L.get(1).toString(),"B");
	}

@Test
public void testCharSequenceIterator() {
	final Iterator<CharSequence> iter =CharSplitter.COMMA.charSequenceIterator(",AA,BBB,,DDDD,,");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next().toString(),"");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next().toString(),"AA");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next().toString(),"BBB");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next().toString(),"");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next().toString(),"DDDD");
	Assert.assertFalse(iter.hasNext());
	}

@Test
public void testStringIterator() {
	final Iterator<String> iter =CharSplitter.COMMA.iterator(",AA,BBB,,DDDD,,");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next(),"");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next(),"AA");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next(),"BBB");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next(),"");
	Assert.assertTrue(iter.hasNext());
	Assert.assertEquals(iter.next(),"DDDD");
	Assert.assertFalse(iter.hasNext());
	}


@Test
public void testStringList01() {
	final List<String> L=CharSplitter.COMMA.splitAsStringList(",AA,BBB,,DDDD,,");
	Assert.assertEquals(L.size(),5);
	Assert.assertEquals(L.get(0),"");
	Assert.assertEquals(L.get(1),"AA");
	Assert.assertEquals(L.get(2),"BBB");
	Assert.assertEquals(L.get(3),"");
	Assert.assertEquals(L.get(4),"DDDD");
	}
@Test
public void testStringList02() {
	final List<String> L=CharSplitter.COMMA.splitAsStringList(",AA,BBB,,DDDD,,",3);
	Assert.assertEquals(L.size(),3);
	Assert.assertEquals(L.get(0),"");
	Assert.assertEquals(L.get(1),"AA");
	Assert.assertEquals(L.get(2),"BBB,,DDDD");
	}
@Test
public void testStringList03() {
	final List<String> L=CharSplitter.COMMA.splitAsStringList("A,B",1000);
	Assert.assertEquals(L.size(),2);
	Assert.assertEquals(L.get(0),"A");
	Assert.assertEquals(L.get(1),"B");
	}

@Test
public void testStringStream() {
	final List<String> L=CharSplitter.COMMA.stream(",AA,BBB,,DDDD,,").collect(Collectors.toList());
	Assert.assertEquals(L.size(),5);
	Assert.assertEquals(L.get(0),"");
	Assert.assertEquals(L.get(1),"AA");
	Assert.assertEquals(L.get(2),"BBB");
	Assert.assertEquals(L.get(3),"");
	Assert.assertEquals(L.get(4),"DDDD");
	}


@Test
public void testStringArray01() {
	final String L[]=CharSplitter.COMMA.split(",AA,BBB,,DDDD,,");
	Assert.assertEquals(L.length,5);
	Assert.assertEquals(L[0],"");
	Assert.assertEquals(L[1],"AA");
	Assert.assertEquals(L[2],"BBB");
	Assert.assertEquals(L[3],"");
	Assert.assertEquals(L[4],"DDDD");
	}
@Test
public void testStringArray02() {
	final String L[]=CharSplitter.COMMA.split(",AA,BBB,,DDDD,,",3);
	Assert.assertEquals(L.length,3);
	Assert.assertEquals(L[0],"");
	Assert.assertEquals(L[1],"AA");
	Assert.assertEquals(L[2],"BBB,,DDDD");
	}
@Test
public void testStringArray03() {
	final String L[]=CharSplitter.COMMA.split("A,B",1000);
	Assert.assertEquals(L.length,2);
	Assert.assertEquals(L[0],"A");
	Assert.assertEquals(L[1],"B");
	}

}
