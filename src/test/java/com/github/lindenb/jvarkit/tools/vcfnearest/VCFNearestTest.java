package com.github.lindenb.jvarkit.tools.vcfnearest;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.variant.infotable.VCFInfoTable;
import com.github.lindenb.jvarkit.variant.infotable.VCFInfoTableModelFactory;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class VCFNearestTest {
	private String createVCF() throws IOException {
		
		String vcf=	"##fileformat=VCFv4.2\n"+
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
			"1\t1000\t.\tG\tC\t0\tPASS\t.\n";
			;
			
		return vcf;
		}
	@Test
	public void testBed() throws Exception {
		final TestSupport support = new TestSupport();
		try {
			String tag = "NEAR";
			Path vin = support.createTmpPath(".vcf");
			Path vout = support.createTmpPath(".vcf");
			Path bed = support.createTmpPath(".bed");
			Files.writeString(vin, createVCF());
			Files.writeString(bed, 
					"1\t999\t1000\tG1\t.\tT1\n"+
					"1\t1000\t1001\tG2\t.\tT2\n"+
					"1\t1001\t1002\tG3\t.\tT3\n"+
					"1\t2000\t3001\tG4\t.\tT4\n"
					);
			
			
			Assert.assertEquals(0,new VCFNearest().instanceMain(new String[] {
					"--tag",tag,
					"--bed",bed.toString(),
					"-C",",,,GENE,,TYPE",
					"-d","100",
					"-n","-1",
					"-o",vout.toString(),
					vin.toString()
					}));
			
			try(VCFIterator iter=new VCFIteratorBuilder().open(vout)) {
				List<VCFInfoTable> tables = new VCFInfoTableModelFactory().parse(iter.getHeader());
				Assert.assertEquals(tables.size(),1);
				VCFInfoTable table=tables.get(0);
				Assert.assertEquals(table.getTag(),tag);

				
				VCFInfoHeaderLine hdr= iter.getHeader().getInfoHeaderLine(tag);
				Assert.assertNotNull(hdr);
				Assert.assertEquals(hdr.getCountType(),VCFHeaderLineCount.UNBOUNDED);
				Assert.assertEquals(hdr.getType(),VCFHeaderLineType.String);
				Assert.assertTrue(iter.hasNext());
				VariantContext ctx = iter.next();
				Assert.assertTrue(ctx.hasAttribute(tag));
				Assert.assertEquals(
						new HashSet<>(ctx.getAttributeAsStringList(tag, "")),
						new HashSet<>(Arrays.asList("G1|T1|0|0","G2|T2|1|2","G3|T3|1|3"))
						);
				}
			}
		catch(Throwable err) {
			Assert.fail("testBed", err);
			}
		finally {
			support.removeTmpFiles();
			}
		}
	
	
	

	}