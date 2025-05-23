package com.github.lindenb.jvarkit.variant.infotable;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VCFInfoTableModelFactoryTest {
    @DataProvider(name = "infoHeaderLines")
    public Object[][] infoHeaderLines() {
        // VCFInfoHeaderLine(ID, Number, Type, Description)
        return new Object[][]{
            {
                new VCFInfoHeaderLine("TABLE", 1, VCFHeaderLineType.String,
                        "This is a table jvarkit.info.bean=(A|B|C)"),
                Arrays.asList("A", "B", "C"),
                "This is a table"
            },
            {
                new VCFInfoHeaderLine("TAB", 1, VCFHeaderLineType.String,
                        "desc jvarkit.info.bean=(X|Y) more"),
                Arrays.asList("X", "Y"),
                "desc more"
            },
            {
                new VCFInfoHeaderLine("NO_TABLE", 1, VCFHeaderLineType.String,
                        "Just a description"),
                null,
                null
            },
            {
                new VCFInfoHeaderLine("WRONG_TYPE", 1, VCFHeaderLineType.Float,
                        "desc jvarkit.info.bean=(A|B)"),
                null,
                null
            }
        };
    }

    @Test(dataProvider = "infoHeaderLines")
    public void testParseInfoHeaderLine(VCFInfoHeaderLine line, List<String> expectedCols, String expectedDesc) {
        VCFInfoTableModelFactory factory = new VCFInfoTableModelFactory();
        Optional<VCFInfoTable> opt = factory.parseInfoHeaderLine(line);
        if (expectedCols == null) {
            Assert.assertFalse(opt.isPresent(), "Expected no table for " + line.getID());
        } else {
            Assert.assertTrue(opt.isPresent(), "Expected table for " + line.getID());
            VCFInfoTable table = opt.get();
            Assert.assertEquals(table.getColumnNames(), expectedCols, "Columns mismatch for " + line.getID());
            Assert.assertTrue(table.getDescription().startsWith(expectedDesc), "Description mismatch for " + line.getID());
        }
    }

    @Test
    public void testEncode() {
        List<String> cols = Arrays.asList("A", "B", "C");
        String encoded = VCFInfoTableModelFactory.encode(cols);
        Assert.assertEquals(encoded, "jvarkit.info.bean=(A|B|C)");
    }
    

    @Test
    public void testParseHeader() {
        VCFInfoHeaderLine line1 = new VCFInfoHeaderLine("T1", 1, VCFHeaderLineType.String,
                "table1 jvarkit.info.bean=(X|Y)");
        VCFInfoHeaderLine line2 = new VCFInfoHeaderLine("T2", 1, VCFHeaderLineType.String,
                "table2 jvarkit.info.bean=(Z)");
        Set<VCFHeaderLine> headerLines = new HashSet<>();
        headerLines.add(line1);
        headerLines.add(line2);
        VCFHeader header = new VCFHeader(headerLines);

        VCFInfoTableModelFactory factory = new VCFInfoTableModelFactory();
        List<VCFInfoTable> tables = factory.parse(header);
        Assert.assertEquals(tables.size(), 2);
        Set<String> tags = new HashSet<>();
        for (VCFInfoTable t : tables) tags.add(t.getTag());
        Assert.assertTrue(tags.contains("T1"));
        Assert.assertTrue(tags.contains("T2"));
    }

    @Test
    public void testNoInfoAttribute() {
        VCFInfoHeaderLine line = new VCFInfoHeaderLine("TABLE", 1, VCFHeaderLineType.String,
                "desc jvarkit.info.bean=(A|B)");
        VCFInfoTableModelFactory factory = new VCFInfoTableModelFactory();
        Optional<VCFInfoTable> opt = factory.parseInfoHeaderLine(line);
        Assert.assertTrue(opt.isPresent());
        VCFInfoTable table = opt.get();
        VariantContext ctx = new VariantContextBuilder("src", "1", 1, 1, Arrays.asList(Allele.create("N", true))).make();
        List<List<String>> rows = table.parse(ctx);
        Assert.assertEquals(rows.size(), 0);
    }
}
