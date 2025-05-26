package com.github.lindenb.jvarkit.bed;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.util.Interval;



public class BedLineTest {

    @Test
    public void testConstructorAndGetters() {
        String[] tokens = {"chr1", "100", "200", "foo"};
        BedLine bedLine = new BedLine(tokens);
        Assert.assertEquals(bedLine.getContig(), "chr1");
        Assert.assertEquals(bedLine.getStart(), 101); // BED is 0-based, Feature is 1-based
        Assert.assertEquals(bedLine.getEnd(), 200);
        Assert.assertEquals(bedLine.getBedStart(), 100);
        Assert.assertEquals(bedLine.getBedEnd(), 200);
        Assert.assertEquals(bedLine.get(0), "chr1");
        Assert.assertEquals(bedLine.get(3), "foo");
        Assert.assertNull(bedLine.get(4));
        Assert.assertEquals(bedLine.getOrDefault(4, "bar"), "bar");
        Assert.assertEquals(bedLine.getColumnCount(), 4);
        Assert.assertEquals(bedLine.join(), "chr1\t100\t200\tfoo");
        Assert.assertEquals(bedLine.join(","), "chr1,100,200,foo");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructor_EmptyContig() {
        String[] tokens = {"", "100", "200"};
        new BedLine(tokens);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructor_BadStart() {
        String[] tokens = {"chr1", "badstart", "200"};
        new BedLine(tokens);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testConstructor_BadEnd() {
        String[] tokens = {"chr1", "100", "badend"};
        new BedLine(tokens);
    }

    @Test
    public void testToInterval() {
        BedLine bedLine = new BedLine(new String[]{"chr2", "9", "20"});
        Interval interval = bedLine.toInterval();
        Assert.assertEquals(interval.getContig(), "chr2");
        Assert.assertEquals(interval.getStart(), 10);
        Assert.assertEquals(interval.getEnd(), 20);
    }

    @Test
    public void testToQueryInterval() {
        BedLine bedLine = new BedLine(new String[]{"chr3", "0", "9"});
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr3", 1000));
        QueryInterval qInt = bedLine.toQueryInterval(dict);
        Assert.assertEquals(qInt.referenceIndex, 0);
        Assert.assertEquals(qInt.start, 1);
        Assert.assertEquals(qInt.end, 9);
    }

    @Test(expectedExceptions = com.github.lindenb.jvarkit.lang.JvarkitException.ContigNotFoundInDictionary.class)
    public void testToQueryInterval_ContigNotFound() {
        BedLine bedLine = new BedLine(new String[]{"chrX", "0", "10"});
        SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", 1000));
        bedLine.toQueryInterval(dict);
    }

    @Test
    public void testRenameContig() {
        String[] tokens = {"chr1", "100", "200"};
        BedLine bedLine = new BedLine(tokens);
        BedLine renamed = bedLine.renameContig("chr2");
        Assert.assertEquals(renamed.getContig(), "chr2");
        Assert.assertEquals(renamed.getStart(), bedLine.getStart());
        Assert.assertEquals(renamed.getEnd(), bedLine.getEnd());
        // Should return same object if contig is unchanged
        Assert.assertSame(bedLine.renameContig("chr1"), bedLine);
    }

    @Test
    public void testIsBedHeader() {
        Assert.assertTrue(BedLine.isBedHeader("# A comment"));
        Assert.assertTrue(BedLine.isBedHeader("track type=bed"));
        Assert.assertTrue(BedLine.isBedHeader("browser position=chr1:1-1000"));
        Assert.assertFalse(BedLine.isBedHeader("chr1\t100\t200"));
    }

    @Test
    public void testToStringArrayAndEqualsHashCode() {
        String[] tokens = {"chr4", "10", "20", "foo"};
        BedLine a = new BedLine(tokens);
        BedLine b = new BedLine(tokens.clone());
        Assert.assertEquals(a, b);
        Assert.assertEquals(a.hashCode(), b.hashCode());
        Assert.assertNotSame(a.toStringArray(), b.toStringArray());
        Assert.assertEquals(a.toStringArray(), b.toStringArray());
        Assert.assertEquals(a.toString(), "chr4\t10\t20\tfoo");
    }

    @Test
    public void testEquals_SelfAndNull() {
        BedLine a = new BedLine(new String[]{"chr1", "100", "101"});
        Assert.assertTrue(a.equals(a));
        Assert.assertFalse(a.equals(null));
        Assert.assertFalse(a.equals(new Object()));
    }

    @DataProvider
    public Object[][] bedHeaderProvider() {
        return new Object[][]{
            {"#foo", true},
            {"track name=foo", true},
            {"browser position=chr1:100-200", true},
            {"chr1\t100\t200", false},
        };
    }

    @Test(dataProvider = "bedHeaderProvider")
    public void testIsBedHeader_DataProvider(String line, boolean expected) {
        Assert.assertEquals(BedLine.isBedHeader(line), expected);
    }
}