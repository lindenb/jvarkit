package com.github.lindenb.jvarkit.variant;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextCodecTest {
	@Test
	public void testLoop() throws Exception {
		final VCFHeader hdr= new VCFHeader();
		final  com.github.lindenb.jvarkit.variant.VariantContextCodec codec=new VariantContextCodec(hdr);
		final List<VariantContext> L=IntStream.range(0, 100)
				.mapToObj(i->new VariantContextBuilder(
					"me",
					"chr"+i, i, i,
					Arrays.asList(Allele.REF_A,Allele.ALT_C)
					).make())
				.collect(Collectors.toList());
		try(ByteArrayOutputStream baos= new ByteArrayOutputStream()) {
			final DataOutputStream dos=new DataOutputStream(baos);
			for(VariantContext ctx: L) codec.encode(dos,ctx);
			baos.flush();
			
			try(DataInputStream dis=new DataInputStream(new ByteArrayInputStream(baos.toByteArray()))) {
				for(int i=0;i< L.size();i++) {
					final VariantContext c1  = L.get(i);
					final VariantContext c2 = codec.decode(dis);
					Assert.assertNotNull(c2);
					Assert.assertEquals(c1.getContig(), c2.getContig());
					Assert.assertEquals(c1.getStart(), c2.getStart());
					Assert.assertEquals(c1.getEnd(), c2.getEnd());
					Assert.assertEquals(c1.getAlleles(), c2.getAlleles());
					}
			}
		}
		
	}
}