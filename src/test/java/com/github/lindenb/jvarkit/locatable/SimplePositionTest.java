package com.github.lindenb.jvarkit.locatable;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Locatable;

public class SimplePositionTest {

    @DataProvider
    public Object[][] validStringProvider() {
        return new Object[][] {
            {"chr1:123", "chr1", 123},
            {"chr2:1", "chr2", 1},
            {"chrX:999", "chrX", 999},
            {"  contigA  :42", "contigA", 42} // test whitespace trimming
        };
    }

    @Test(dataProvider = "validStringProvider")
    public void testConstructorFromString(String input, String expectedContig, int expectedPos) {
        SimplePosition pos = new SimplePosition(input);
        Assert.assertEquals(pos.getContig(), expectedContig);
        Assert.assertEquals(pos.getPosition(), expectedPos);
        Assert.assertEquals(pos.getStart(), expectedPos);
        Assert.assertEquals(pos.getEnd(), expectedPos);
    }

    @DataProvider
    public Object[][] invalidStringProvider() {
        return new Object[][] {
            {null},
            {""},
            {"chr1"},
            {"chr1:"},
            {":123"},
            {"   :456"},
            {"chr1:abc"}
        };
    }

    @Test(dataProvider = "invalidStringProvider", expectedExceptions = IllegalArgumentException.class)
    public void testConstructorFromStringBadCases(String input) {
        new SimplePosition(input);
    }

    @Test
    public void testConstructorFromContigAndPos() {
        SimplePosition pos = new SimplePosition("chr3", 888);
        Assert.assertEquals(pos.getContig(), "chr3");
        Assert.assertEquals(pos.getPosition(), 888);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructorNullContig() {
        new SimplePosition(null, 100);
    }

    @Test
    public void testLengthAndGetLengthOnReference() {
        SimplePosition pos = new SimplePosition("chr1", 123);
        Assert.assertEquals(pos.length(), 1);
        Assert.assertEquals(pos.getLengthOnReference(), 1);
    }

    @Test
    public void testEqualsAndHashCode() {
        SimplePosition p1 = new SimplePosition("chr1", 123);
        SimplePosition p2 = new SimplePosition("chr1", 123);
        SimplePosition p3 = new SimplePosition("chr1", 124);
        SimplePosition p4 = new SimplePosition("chr2", 123);

        Assert.assertEquals(p1, p2);
        Assert.assertEquals(p1.hashCode(), p2.hashCode());
        Assert.assertNotEquals(p1, p3);
        Assert.assertNotEquals(p1, p4);
        Assert.assertNotEquals(p1, null);
        Assert.assertNotEquals(p1, "chr1:123");
    }

    @Test
    public void testCompareTo() {
        SimplePosition p1 = new SimplePosition("chr1", 123);
        SimplePosition p2 = new SimplePosition("chr1", 124);
        SimplePosition p3 = new SimplePosition("chr2", 100);

        Assert.assertTrue(p1.compareTo(p2) < 0);
        Assert.assertTrue(p2.compareTo(p1) > 0);
        Assert.assertTrue(p1.compareTo(p3) < 0);
        Assert.assertTrue(p3.compareTo(p1) > 0);
        Assert.assertEquals(p1.compareTo(new SimplePosition("chr1", 123)), 0);
    }

    @Test
    public void testToString() {
        SimplePosition pos = new SimplePosition("chrM", 7);
        Assert.assertEquals(pos.toString(), "chrM:7");
    }

    @Test
    public void testRenameContig() {
        SimplePosition pos = new SimplePosition("foo", 3);
        SimplePosition renamed = pos.renameContig("bar");
        Assert.assertEquals(renamed.getContig(), "bar");
        Assert.assertEquals(renamed.getPosition(), 3);
        Assert.assertNotSame(renamed, pos);

        // If same contig, should return the same object
        Assert.assertSame(pos.renameContig("foo"), pos);
    }

    @Test
    public void testExtendPositive() {
        SimplePosition pos = new SimplePosition("chr1", 10);
        Locatable interval = pos.extend(5);
        Assert.assertTrue(interval instanceof SimpleInterval);
        Assert.assertEquals(interval.getContig(), "chr1");
        Assert.assertEquals(interval.getStart(), 5); // max(1, 10-5)
        Assert.assertEquals(interval.getEnd(), 15);
    }

    @Test
    public void testExtendZero() {
        SimplePosition pos = new SimplePosition("chr1", 10);
        Locatable result = pos.extend(0);
        Assert.assertSame(result, pos);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testExtendNegative() {
        SimplePosition pos = new SimplePosition("chr1", 10);
        pos.extend(-2);
    }
}