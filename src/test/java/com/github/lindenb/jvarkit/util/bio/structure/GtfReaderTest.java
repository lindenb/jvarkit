package com.github.lindenb.jvarkit.util.bio.structure;

import java.io.IOException;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;

@AlsoTest(PeptideSequenceTest.class)
public class GtfReaderTest {
private final TestSupport support =new TestSupport();	
	
@Test
void test01() throws IOException {
	try {
		String path = support.resource("Homo_sapiens.GRCh37.87.gtf.gz");
		GtfReader gf = new GtfReader(path);
		final List<Gene> genes = gf.getAllGenes();
		gf.close();
		Assert.assertEquals(genes.size(), 3);
		
		Gene g  = genes.stream().filter(G->G.getId().equals("ENSG00000134250")).findFirst().orElseThrow(()->new IOException("Cannot get ENSG00000134250 "));
		Assert.assertEquals(g.getStart(), 120454176);
		Assert.assertEquals(g.getEnd(), 120612240);
		Assert.assertEquals(g.getContig(), "1");
		Assert.assertEquals(g.getStrand(),'-');
		Assert.assertTrue(g.isNegativeStrand());
		Assert.assertFalse(g.isPositiveStrand());
		Assert.assertEquals(g.getTranscripts().size(),6);
		
		Transcript tr  = g.getTranscripts().stream().filter(G->G.getId().equals("ENST00000602566")).findFirst().orElseThrow(()->new IOException("Cannot get ENST00000602566 "));
		Assert.assertEquals(tr.getStart(), 120535404);
		Assert.assertEquals(tr.getEnd(), 120596839);
		Assert.assertEquals(tr.getTxStart(), tr.getStart());
		Assert.assertEquals(tr.getTxEnd(), tr.getEnd());
		
		}
	finally
		{
		
		}
	}
}
