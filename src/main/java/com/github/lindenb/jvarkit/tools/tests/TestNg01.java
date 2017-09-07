/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


*/
package com.github.lindenb.jvarkit.tools.tests;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.testng.Assert;
import org.testng.annotations.*;

import com.github.lindenb.jvarkit.tools.bam2graphics.Bam2Raster;
import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk;
import com.github.lindenb.jvarkit.tools.biostar.Biostar86480;
import com.github.lindenb.jvarkit.tools.burden.VcfInjectPedigree;
import com.github.lindenb.jvarkit.tools.burden.VcfLoopOverGenes;
import com.github.lindenb.jvarkit.tools.groupbygene.GroupByGene;
import com.github.lindenb.jvarkit.tools.misc.VCFPolyX;
import com.github.lindenb.jvarkit.tools.misc.VcfCreateDictionary;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.misc.VcfSetSequenceDictionary;
import com.github.lindenb.jvarkit.tools.misc.VcfTail;
import com.github.lindenb.jvarkit.tools.misc.VcfToHilbert;
import com.github.lindenb.jvarkit.tools.misc.VcfToSvg;
import com.github.lindenb.jvarkit.tools.misc.VcfToTable;
import com.github.lindenb.jvarkit.tools.sam2tsv.Sam2Tsv;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcffixindels.VCFFixIndels;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

class TestNg01 {
	static final File TEST_RESULTS_DIR= new File("test-results");
	static final File JETER_VCF = new File(TEST_RESULTS_DIR,"jeter.vcf");
	static final String TOY_FA="src/test/resources/toy.fa";
	static final String TOY_VCF_GZ="src/test/resources/toy.vcf.gz";
	static final String TOY_BAM="src/test/resources/toy.bam";
	static final String TOY_DICT="src/test/resources/toy.dict";
	static final String VCF01 = "src/test/resources/test_vcf01.vcf";
	static final String PED01 = "src/test/resources/test_vcf01.ped";
	static final String KNOWN_GENES01 = "src/test/resources/test_vcf01.knownGenes.txt.gz";
	@BeforeClass
    public void setup() throws IOException {
		TEST_RESULTS_DIR.mkdirs();
		}
	
	@DataProvider(name = "all_vcfs")
	public Object[][]  all_vcfs() {
		return new Object[][] {{TOY_VCF_GZ},{VCF01}}	;
		}
	
	 static Stream<VariantContext> streamVcf(final File f) {
		final VCFFileReader r = new VCFFileReader(f,false);
		final CloseableIterator<VariantContext> iter = r.iterator();
		return StreamSupport.stream(new IterableAdapter<VariantContext>(iter).spliterator(), false).onClose(()->{iter.close();r.close();});
		}
	 static Stream<VariantContext> streamJeterVcf() {
		return streamVcf(JETER_VCF);
		}
	
	
    @Test(dataProvider="all_vcfs")
    public void testVcf2Table(final String vcfPath) {
    	File output = new File(TEST_RESULTS_DIR,"jeter.txt");
        Assert.assertEquals(0,new VcfToTable().instanceMain(new String[]{
        		"-o",output.getPath(),
        		vcfPath
        	}));
        Assert.assertTrue( output.exists());
    	}
    @Test
    public void testVcfSetSequenceDictionary() throws IOException {
    	File output = new File(TEST_RESULTS_DIR,"jeter.vcf");
    	File dict = new File(TEST_RESULTS_DIR,"jeter.dict");
    	PrintWriter pw = new PrintWriter(dict);
    	BufferedReader br=new BufferedReader(new FileReader(TOY_DICT));
    	br.lines().map(L->L.replace("\tSN:", "\tSN:chr")).forEach(L->pw.println(L));
    	br.close();
    	pw.flush();pw.close();
    	
        Assert.assertEquals(0,new VcfSetSequenceDictionary().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-R",dict.getPath(),
        		TOY_VCF_GZ
        	}));
        Assert.assertTrue( output.exists());
    	}
    @Test(dataProvider="all_vcfs")
    public void testVcfCreateSequenceDictionary(final String vcfPath) throws IOException{
    	File output = new File(TEST_RESULTS_DIR,"jeter.dict");
    
        Assert.assertEquals(0,new VcfCreateDictionary().instanceMain(new String[]{
        		"-o",output.getPath(),
        		vcfPath
        	}));
        Assert.assertTrue( output.exists());
    	}
    @Test
    public void testSam2Tsv() throws IOException{
    	File output = new File(TEST_RESULTS_DIR,"jeter.txt");
    
        Assert.assertEquals(0,new Sam2Tsv().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-R",TOY_FA,
        		TOY_BAM
        	}));
        Assert.assertTrue( output.exists());
    	}
    @Test
    public void testGroupByGene() throws IOException{
    	File output = new File(TEST_RESULTS_DIR,"jeter.txt");
    
        Assert.assertEquals(0,new GroupByGene().instanceMain(new String[]{
        		"-o",output.getPath(),
        		VCF01
        	}));
        Assert.assertTrue( output.exists());
    	}

    @Test(dataProvider="all_vcfs")
    public void testBioAlcidaeJdkVcf(final String vcfPath) throws IOException{
    	File output = new File(TEST_RESULTS_DIR,"jeter.txt");
    
        Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-e","print(stream().count());",
        		vcfPath
        	}));
        Assert.assertTrue( output.exists());
    	}

    
    @Test(dataProvider="all_vcfs")
    public void testVcfFiltrerJdkVcf(final String vcfPath) throws IOException{
        Assert.assertEquals(0,new VcfFilterJdk().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-e","return variant.getStart()%2==0;",
        		vcfPath
        	}));
        Assert.assertTrue( JETER_VCF.exists());
        Assert.assertEquals(streamVcf(JETER_VCF).filter(V->V.getStart()%2!=0).count(),0L);
    	}

    public void testSamJdk() throws IOException{
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.sam");
    	Assert.assertEquals(0,new VcfFilterJdk().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-e","return record.getStart()%2==0;",
        		TOY_BAM
        	}));
        Assert.assertTrue( tmp.exists());
    	}

    public void testBamToRaster() throws IOException{
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.png");
    	Assert.assertEquals(0,new Bam2Raster().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-R",TOY_FA,
        		"-r","seq:1-100",
        		TOY_BAM
        	}));
        Assert.assertTrue( tmp.exists());
    	}
    
    @Test
    public void testVcfPolyX() throws IOException{
        Assert.assertEquals(0,new VCFPolyX().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-R",TOY_FA,
        		"-n","1",
        		TOY_VCF_GZ
        	}));
        Assert.assertTrue( JETER_VCF.exists());
    	}
    @Test(dataProvider="all_vcfs")
    public void testVcfHead(final String vcfPath) throws IOException{    
        Assert.assertEquals(0,new VcfHead().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-n","1",
        		TOY_VCF_GZ
        	}));
        Assert.assertTrue( JETER_VCF.exists());
        Assert.assertEquals(streamVcf(JETER_VCF).count(),1L);
    	}
    @Test(dataProvider="all_vcfs")
    public void testVcfTail(final String vcfPath) throws IOException{    
        Assert.assertEquals(0,new VcfTail().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-n","1",
        		vcfPath
        	}));
        Assert.assertEquals(streamVcf(JETER_VCF).count(),1L);
    	}
    @Test
    public void testVcfInjectPed() throws IOException{    
        Assert.assertEquals(0,new VcfInjectPedigree().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"--pedigree",PED01,
        		VCF01
        	}));
    	}
    @Test
    public void testVcfToSvg() throws IOException{    
    	File svg = new File(TEST_RESULTS_DIR,"jeter.__SEGMENT__.svg");
    	Assert.assertTrue(svg.getName().contains("__SEGMENT__"));
    	Assert.assertEquals(0,new VcfToSvg().instanceMain(new String[]{
        		"-o",svg.getPath(),
        		"--stopAfterFirst",
        		"-k",KNOWN_GENES01,
        		VCF01}));
    	// Assert.assertTrue( svg.exists()); NO, name is changed
    	}
    @Test(dataProvider="all_vcfs")
    public void testVcfToHilbert(final String vcfPath) throws IOException{    
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.png");
    	Assert.assertEquals(0,new VcfToHilbert().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		vcfPath}));
    	Assert.assertTrue( tmp.exists());
    	}
    
    
    @Test(dataProvider="all_vcfs")
    public void testVcfFixIndels(final String vcfPath) throws IOException{    
    	Assert.assertEquals(0,new VCFFixIndels().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		vcfPath}));
    	}

    
    @Test
    public void testVcfLoopOverGenes() throws IOException{    
    	final String prefix="vlog";
    	File genes = new File(TEST_RESULTS_DIR,"jeter.genes.txt");
    	Assert.assertEquals(0,new VcfLoopOverGenes().instanceMain(new String[]{
        		"-o",genes.getPath(),
        		"-p",prefix,
        		"--snpEffNoIntergenic",
        		VCF01
        		}));
    	Assert.assertTrue( genes.exists());
    	
    	Assert.assertEquals(0,new VcfLoopOverGenes().instanceMain(new String[]{
        		"-o",TEST_RESULTS_DIR.getPath(),
        		"-g",genes.getPath(),
        		VCF01
        		}));

    	
    	Assert.assertEquals(0,new VcfLoopOverGenes().instanceMain(new String[]{
        		"-o",genes.getPath(),
        		"--splitMethod","VariantSlidingWindow",
        		"--variantsWinCount","10",
        		"--variantsWinShift","10",
        		"-p",prefix,
        		"--snpEffNoIntergenic",
        		VCF01
        		}));
    	
    	Assert.assertEquals(0,new VcfLoopOverGenes().instanceMain(new String[]{
        		"-o",TEST_RESULTS_DIR.getPath(),
        		"-g",genes.getPath(),
        		VCF01
        		}));
    	Arrays.asList(TEST_RESULTS_DIR.listFiles()).stream().
    		filter(F->F.getName().startsWith(prefix) && F.getName().endsWith(".vcf")).
    		forEach(F->F.delete());
    	}
    @Test
    public void testBiostar86480() throws IOException{    
    	File tmp = new File(TEST_RESULTS_DIR, "tmp.txt");
    	Assert.assertEquals(0,new Biostar86480().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-E","EcoRI",
        		"-E","BamHI",
        		TOY_FA
        		}));
    	}
    @Test
    public void testVcfFilterSo() throws IOException{   
    	final AnnPredictionParser parser = new AnnPredictionParserFactory().createDefaultParser();
    	final SequenceOntologyTree tree = SequenceOntologyTree.getInstance();
    	String acn = "SO:0001583";
    	final SequenceOntologyTree.Term term =tree.getTermByAcn(acn);
    	final Set<SequenceOntologyTree.Term> terms = term.getAllDescendants();
    	Assert.assertNotNull(term);
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-A",acn,
        		VCF01
        		}));
    	streamVcf(JETER_VCF).forEach(V->{
    		Assert.assertTrue(parser.getPredictions(V).stream().
    			flatMap(P->P.getSOTerms().stream()).
				anyMatch(T->terms.contains(V)));
    	});
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-A",acn,
        		"--invert",
        		VCF01
        		}));
    	streamVcf(JETER_VCF).forEach(V->{
    		Assert.assertFalse(parser.getPredictions(V).stream().
    			flatMap(P->P.getSOTerms().stream()).
				anyMatch(T->terms.contains(V)));
    	});

    	
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-A",acn,
        		"--rmatt",
        		VCF01
        		}));
    	Assert.assertTrue(streamVcf(JETER_VCF).findAny().isPresent());
    	}

	}
