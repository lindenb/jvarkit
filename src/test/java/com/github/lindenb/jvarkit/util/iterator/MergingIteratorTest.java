package com.github.lindenb.jvarkit.util.iterator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class MergingIteratorTest {
@Test
public void test1() {
	final MergingIterator<Integer> iter = new MergingIterator<>(
			Integer::compare,
			Arrays.asList(
					Arrays.asList(1,3,5,7,9).iterator(),
					Arrays.asList(0,2,4,6,8,10).iterator(),
					Arrays.asList(1,3,5,7,9).iterator(),
					Arrays.asList(0,2,4,6,8,10).iterator(),
					Arrays.asList(11,11,12,12).iterator()
			));
	for(int i=0;i<=12;i++)
		{
		for(int n=0;n<2;++n) {
			Assert.assertTrue(iter.hasNext());
			int j= iter.next();
			Assert.assertEquals(i,j);
			}
		}
	Assert.assertFalse(iter.hasNext());
	iter.close();
	}
@Test(expectedExceptions= {IllegalStateException.class})
public void test2() {
	final MergingIterator<Integer> iter = new MergingIterator<>(
			Integer::compare,
			Arrays.asList(
					Arrays.asList(10,9).iterator(),
					Arrays.asList(4,2).iterator()
			));
	while(iter.hasNext()) iter.next();
	iter.close();
	}

@Test
public void test3() {
	SAMSequenceDictionary dict= new SAMSequenceDictionary();
	dict.addSequence(new SAMSequenceRecord("1",10000));
	SAMFileHeader header = new SAMFileHeader();
	header.setSequenceDictionary(dict);
	final List<SAMRecord> L = new ArrayList<>();
	DefaultSAMRecordFactory srf = new DefaultSAMRecordFactory();
	SAMRecord r1= srf.createSAMRecord(header);
	
	r1.setReadName("R1");
	r1.setReferenceName("1");
	r1.setAlignmentStart(1);
	r1.setFlags(99);
	r1.setReadString("GAATTC");
	r1.setBaseQualityString("222222");
	r1.setCigarString("6M");
	L.add(r1);
	
	r1= srf.createSAMRecord(header);
	
	r1.setReadName("R2");
	r1.setReferenceName("1");
	r1.setAlignmentStart(83);
	r1.setFlags(83);
	r1.setReadString("GAATTC");
	r1.setBaseQualityString("222222");
	r1.setCigarString("6M");
	L.add(r1);
	
	final MergingIterator<SAMRecord> iter = new MergingIterator<>(
			new SAMRecordCoordinateComparator(),
			Collections.singletonList(L.iterator())
			);
	while(iter.hasNext()) iter.next();
	iter.close();
	}

}
