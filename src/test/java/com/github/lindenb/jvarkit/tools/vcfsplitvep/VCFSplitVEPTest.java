package com.github.lindenb.jvarkit.tools.vcfsplitvep;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashSet;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.variant.vcf.BcfIteratorBuilder;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

public class VCFSplitVEPTest {
	private String createVCF(String csq,String...infos) throws IOException {
		String info="CSQ="+Arrays.stream(infos).map(S->"a|b|"+S+"|d").collect(Collectors.joining(","));
		
		String vcf=	"##fileformat=VCFv4.2\n"+
			"##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: x1|x2|"+csq+"|x4\">\n"+
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
			"1\t1\t.\tG\tC\t0\tPASS\t"+info+"\n";
			;
			
		return vcf;
		}
	@Test
	public void testMin() throws Exception {
		Path vin = null;
		Path vout = null;
		try {
			vin= Files.createTempFile("tmp", ".vcf");
			vout= Files.createTempFile("tmp", ".vcf");
			Files.writeString(vin, createVCF("maf","1","10","12","3","","1.2"));
			Assert.assertEquals(0,new VCFSplitVEP().instanceMain(new String[] {"--tags","maf:Float:min","-o",vout.toString(),vin.toString()}));
			String id = "VEP_maf_min";
			try(VCFIterator iter=new BcfIteratorBuilder().open(vout)) {
				VCFInfoHeaderLine hdr= iter.getHeader().getInfoHeaderLine(id);
				Assert.assertNotNull(hdr);
				Assert.assertEquals(hdr.getCount(),1);
				Assert.assertEquals(hdr.getType(),VCFHeaderLineType.Float);
				Assert.assertTrue(iter.hasNext());
				VariantContext ctx = iter.next();
				Assert.assertTrue(ctx.hasAttribute(id));
				Assert.assertEquals(ctx.getAttributeAsDouble(id,-1.0),1.0);
				}
			}
		catch(Throwable err) {
			Assert.fail("min", err);
			}
		finally {
			if(vin!=null) IOUtil.deletePaths(vin);
			if(vout!=null) IOUtil.deletePaths(vout);
			}
		}
	
	@Test
	public void testMax() throws Exception {
		Path vin = null;
		Path vout = null;
		try {
			vin= Files.createTempFile("tmp", ".vcf");
			vout= Files.createTempFile("tmp", ".vcf");
			Files.writeString(vin, createVCF("maf","1","100","12","3","","1"));
			Assert.assertEquals(0,new VCFSplitVEP().instanceMain(new String[] {"--tags","maf:Integer:max","-o",vout.toString(),vin.toString()}));
			String id = "VEP_maf_max";
			try(VCFIterator iter=new BcfIteratorBuilder().open(vout)) {
				VCFInfoHeaderLine hdr= iter.getHeader().getInfoHeaderLine(id);
				Assert.assertNotNull(hdr);
				Assert.assertEquals(hdr.getCount(),1);
				Assert.assertEquals(hdr.getType(),VCFHeaderLineType.Integer);
				Assert.assertTrue(iter.hasNext());
				VariantContext ctx = iter.next();
				Assert.assertTrue(ctx.hasAttribute(id));
				Assert.assertEquals(ctx.getAttributeAsInt(id,-1),100);
				}
			}
		catch(Throwable err) {
			Assert.fail("max", err);
			}
		finally {
			if(vin!=null) IOUtil.deletePaths(vin);
			if(vout!=null) IOUtil.deletePaths(vout);
			}
		}
	
	
	@Test
	public void testRandom() throws Exception {
		Path vin = null;
		Path vout = null;
		try {
			vin= Files.createTempFile("tmp", ".vcf");
			vout= Files.createTempFile("tmp", ".vcf");
			Files.writeString(vin, createVCF("maf","1","100","12","3","","1"));
			Assert.assertEquals(0,new VCFSplitVEP().instanceMain(new String[] {"--tags","maf:Integer:random","-o",vout.toString(),vin.toString()}));
			String id = "VEP_maf_random";
			try(VCFIterator iter=new BcfIteratorBuilder().open(vout)) {
				VCFInfoHeaderLine hdr= iter.getHeader().getInfoHeaderLine(id);
				Assert.assertNotNull(hdr);
				Assert.assertEquals(hdr.getCount(),1);
				Assert.assertEquals(hdr.getType(),VCFHeaderLineType.Integer);
				Assert.assertTrue(iter.hasNext());
				VariantContext ctx = iter.next();
				Assert.assertTrue(ctx.hasAttribute(id));
				Assert.assertTrue(new HashSet<>(Arrays.asList(1,100,12,3,1)).contains(ctx.getAttributeAsInt(id,-1)));
				}
			}
		catch(Throwable err) {
			Assert.fail("random", err);
			}
		finally {
			if(vin!=null) IOUtil.deletePaths(vin);
			if(vout!=null) IOUtil.deletePaths(vout);
			}
		}
	
	@Test
	public void testUniq() throws Exception {
		Path vin = null;
		Path vout = null;
		try {
			vin= Files.createTempFile("tmp", ".vcf");
			vout= Files.createTempFile("tmp", ".vcf");
			Files.writeString(vin, createVCF("maf","A","A","B","C","","D","C"));
			Assert.assertEquals(0,new VCFSplitVEP().instanceMain(new String[] {"--tags","maf:String:uniq","-o",vout.toString(),vin.toString()}));
			String id = "VEP_maf_uniq";
			try(VCFIterator iter=new BcfIteratorBuilder().open(vout)) {
				VCFInfoHeaderLine hdr= iter.getHeader().getInfoHeaderLine(id);
				Assert.assertNotNull(hdr);
				Assert.assertEquals(hdr.getCountType(),VCFHeaderLineCount.UNBOUNDED);
				Assert.assertEquals(hdr.getType(),VCFHeaderLineType.String);
				Assert.assertTrue(iter.hasNext());
				VariantContext ctx = iter.next();
				Assert.assertTrue(ctx.hasAttribute(id));
				Assert.assertEquals(new HashSet<>(ctx.getAttributeAsStringList(id,"")),new HashSet<>(Arrays.asList("A","B","C","D")));
				}
			}
		catch(Throwable err) {
			Assert.fail("uniq", err);
			}
		finally {
			if(vin!=null) IOUtil.deletePaths(vin);
			if(vout!=null) IOUtil.deletePaths(vout);
			}
		}
	}