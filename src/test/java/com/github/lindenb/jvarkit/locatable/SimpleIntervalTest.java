package com.github.lindenb.jvarkit.locatable;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import htsjdk.samtools.util.Locatable;

public class SimpleIntervalTest {

    @DataProvider
    public Object[][] fromStringProvider() {
        return new Object[][]{
            {"chr1:100-200", "chr1", 100, 200},
            {"chr2:150", "chr2", 150, 150},
            {"chr3:200+10", "chr3", 190, 210},
            {"chr4:10-10", "chr4", 10, 10},
            {"chr5:5,000-6,000", "chr5", 5000, 6000}
        };
    }

    @Test(dataProvider="fromStringProvider")
    public void testFromString(String str, String contig, int start, int end) {
        SimpleInterval si = new SimpleInterval(str);
        Assert.assertEquals(si.getContig(), contig);
        Assert.assertEquals(si.getStart(), start);
        Assert.assertEquals(si.getEnd(), end);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNullStringConstructor() {
        new SimpleInterval((String) null);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNoColon() {
        new SimpleInterval("no_colon");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeExtend() {
        SimpleInterval si = new SimpleInterval("chr1:100-200");
        si.extend(-5);
    }

    @Test
    public void testExtend() {
        SimpleInterval si = new SimpleInterval("chr1:100-200");
        SimpleInterval extended = si.extend(10);
        Assert.assertEquals(extended.getStart(), 90);
        Assert.assertEquals(extended.getEnd(), 210);
        Assert.assertEquals(extended.getContig(), "chr1");
        Assert.assertNotEquals(extended, si);
        Assert.assertEquals(si.extend(0), si);
    }

    @Test
    public void testRenameContig() {
        SimpleInterval si = new SimpleInterval("chr1:100-200");
        SimpleInterval renamed = si.renameContig("chrX");
        Assert.assertEquals(renamed.getContig(), "chrX");
        Assert.assertEquals(renamed.getStart(), 100);
        Assert.assertEquals(renamed.getEnd(), 200);
        Assert.assertNotEquals(renamed, si);
        Assert.assertEquals(si.renameContig("chr1"), si);
    }

    @Test
    public void testMergeOverlap() {
        SimpleInterval i1 = new SimpleInterval("chr1:100-150");
        SimpleInterval i2 = new SimpleInterval("chr1:120-200");
        SimpleInterval merged = i1.merge(i2);
        Assert.assertEquals(merged.getContig(), "chr1");
        Assert.assertEquals(merged.getStart(), 100);
        Assert.assertEquals(merged.getEnd(), 200);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testMergeNonOverlap() {
        SimpleInterval i1 = new SimpleInterval("chr1:100-110");
        SimpleInterval i2 = new SimpleInterval("chr1:120-130");
        i1.merge(i2);
    }

    @Test
    public void testMergeCollection() {
        List<SimpleInterval> intervals = Arrays.asList(
                new SimpleInterval("chr1:100-150"),
                new SimpleInterval("chr1:120-180"),
                new SimpleInterval("chr1:200-210"),
                new SimpleInterval("chr1:209-215")
        );
        List<SimpleInterval> merged = SimpleInterval.mergeCollection(intervals);
        Assert.assertEquals(merged.size(), 2);
        Assert.assertEquals(merged.get(0).getStart(), 100);
        Assert.assertEquals(merged.get(0).getEnd(), 180);
        Assert.assertEquals(merged.get(1).getStart(), 200);
        Assert.assertEquals(merged.get(1).getEnd(), 215);
    }

    @Test
    public void testGetCenter() {
        SimpleInterval si = new SimpleInterval("chr1:10-19");
        SimplePosition center = si.getCenter();
        Assert.assertEquals(center.getContig(), "chr1");
        Assert.assertEquals(center.getPosition(), 15);
    }

    @Test
    public void testIntersectionLength() {
        SimpleInterval i1 = new SimpleInterval("chr1:100-200");
        SimpleInterval i2 = new SimpleInterval("chr1:150-250");
        Assert.assertEquals(i1.getIntersectionLength(i2), 51); // 150-200 inclusive = 51 bp
        Assert.assertEquals(i2.getIntersectionLength(i1), 51);
        SimpleInterval i3 = new SimpleInterval("chr2:100-200");
        Assert.assertEquals(i1.getIntersectionLength(i3), 0);
    }

    @Test
    public void testEqualsAndCompareTo() {
        SimpleInterval i1 = new SimpleInterval("chr1:100-200");
        SimpleInterval i2 = new SimpleInterval("chr1:100-200");
        SimpleInterval i3 = new SimpleInterval("chr1:150-200");
        Assert.assertEquals(i1, i2);
        Assert.assertNotEquals(i1, i3);
        Assert.assertTrue(i1.compareTo(i2) == 0);
        Assert.assertTrue(i1.compareTo(i3) < 0);
        Assert.assertTrue(i3.compareTo(i1) > 0);
    }

    @Test
    public void testFromLocatable() {
        Locatable loc = new Locatable() {
            @Override public String getContig() { return "chrZ"; }
            @Override public int getStart() { return 10; }
            @Override public int getEnd() { return 20; }
        };
        SimpleInterval si = new SimpleInterval(loc);
        Assert.assertEquals(si.getContig(), "chrZ");
        Assert.assertEquals(si.getStart(), 10);
        Assert.assertEquals(si.getEnd(), 20);
    }

    @Test
    public void testLength() {
        SimpleInterval si = new SimpleInterval("chr1:100-105");
        Assert.assertEquals(si.length(), 6); // inclusive: 105-100+1
    }

    @Test
    public void testMergeCollectionSingleton() {
        List<SimpleInterval> intervals = Collections.singletonList(new SimpleInterval("chr1:100-200"));
        List<SimpleInterval> merged = SimpleInterval.mergeCollection(intervals);
        Assert.assertEquals(merged.size(), 1);
        Assert.assertEquals(merged.get(0).getStart(), 100);
        Assert.assertEquals(merged.get(0).getEnd(), 200);
    }
}
