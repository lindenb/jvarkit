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
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.testng.Assert;
import org.testng.annotations.*;

import com.github.lindenb.jvarkit.tools.bam2graphics.Bam2Raster;
import com.github.lindenb.jvarkit.tools.bam2graphics.LowResBam2Raster;
import com.github.lindenb.jvarkit.tools.bam2wig.Bam2Wig;
import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidaeJdk;
import com.github.lindenb.jvarkit.tools.biostar.Biostar59647;
import com.github.lindenb.jvarkit.tools.biostar.Biostar86480;
import com.github.lindenb.jvarkit.tools.burden.CaseControlCanvas;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFilterExac;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFilterGenes;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFisherH;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenFisherV;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenMAF;
import com.github.lindenb.jvarkit.tools.burden.VcfBurdenRscriptV;
import com.github.lindenb.jvarkit.tools.burden.VcfFilterNotInPedigree;
import com.github.lindenb.jvarkit.tools.burden.VcfInjectPedigree;
import com.github.lindenb.jvarkit.tools.burden.VcfLoopOverGenes;
import com.github.lindenb.jvarkit.tools.burden.VcfMoveFiltersToInfo;
import com.github.lindenb.jvarkit.tools.calling.MiniCaller;
import com.github.lindenb.jvarkit.tools.fastq.FastqShuffle;
import com.github.lindenb.jvarkit.tools.gnomad.VcfGnomad;
import com.github.lindenb.jvarkit.tools.groupbygene.GroupByGene;
import com.github.lindenb.jvarkit.tools.misc.BamToSql;
import com.github.lindenb.jvarkit.tools.misc.ConvertVcfChromosomes;
import com.github.lindenb.jvarkit.tools.misc.FindAVariation;
import com.github.lindenb.jvarkit.tools.misc.FindAllCoverageAtPosition;
import com.github.lindenb.jvarkit.tools.misc.FixVcfMissingGenotypes;
import com.github.lindenb.jvarkit.tools.misc.Gff2KnownGene;
import com.github.lindenb.jvarkit.tools.misc.PadEmptyFastq;
import com.github.lindenb.jvarkit.tools.misc.VCFPolyX;
import com.github.lindenb.jvarkit.tools.misc.VcfCreateDictionary;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.misc.VcfMultiToOneAllele;
import com.github.lindenb.jvarkit.tools.misc.VcfNoCallToHomRef;
import com.github.lindenb.jvarkit.tools.misc.VcfSetSequenceDictionary;
import com.github.lindenb.jvarkit.tools.misc.VcfTail;
import com.github.lindenb.jvarkit.tools.misc.VcfToHilbert;
import com.github.lindenb.jvarkit.tools.misc.VcfToSvg;
import com.github.lindenb.jvarkit.tools.misc.VcfToTable;
import com.github.lindenb.jvarkit.tools.ngsfiles.NgsFilesSummary;
import com.github.lindenb.jvarkit.tools.onesamplevcf.VcfMultiToOne;
import com.github.lindenb.jvarkit.tools.sam2tsv.Sam2Tsv;
import com.github.lindenb.jvarkit.tools.sam4weblogo.SAM4WebLogo;
import com.github.lindenb.jvarkit.tools.samjs.SamJdk;
import com.github.lindenb.jvarkit.tools.sortvcfonref.SortVcfOnInfo;
import com.github.lindenb.jvarkit.tools.vcf2xml.Vcf2Xml;
import com.github.lindenb.jvarkit.tools.vcfamalgation.VcfXmlAmalgamation;
import com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig;
import com.github.lindenb.jvarkit.tools.vcfcmp.VcfCompareCallers;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcffixindels.VCFFixIndels;
import com.github.lindenb.jvarkit.tools.vcfrebase.VcfRebase;
import com.github.lindenb.jvarkit.tools.vcfstats.VcfStats;
import com.github.lindenb.jvarkit.tools.vcfstripannot.VCFStripAnnotations;
import com.github.lindenb.jvarkit.tools.vcftrios.VCFTrios;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

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
	
	private Properties properties = new Properties();
	
	private static class VcfBuilder
		{
		private final VCFHeader header =new VCFHeader();
		private final List<VariantContext> variants = new ArrayList<>();
		VcfBuilder()
			{
			
			}
		public VcfBuilder add(final VariantContext ctx)
			{
			this.variants.add(ctx);
			return this;
			}
		public VcfBuilder write(final File out)
			{
			final VariantContextWriter w= new VariantContextWriterBuilder().
					setOutputFile(out).
					build();
			w.writeHeader(this.header);
			for(final VariantContext ctx:this.variants) w.add(ctx);
			w.close();
			return this;
			}

		}
	
	@BeforeClass
    public void setup() throws IOException {
		TEST_RESULTS_DIR.mkdirs();
		final File propFile = new File(TEST_RESULTS_DIR,"properties.xml");
		if(propFile.exists()) {
			FileInputStream in = new FileInputStream(propFile);
			this.properties.loadFromXML(in);
			in.close();
			}
		}
	
	@DataProvider(name = "all_vcfs")
	public Object[][]  all_vcfs() {
		return new Object[][] {
			{TOY_VCF_GZ},
			{VCF01},
			{"src/test/resources/ExAC.r1.sites.vep.vcf.gz"},
			{"src/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz"},
			{"src/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz"}
			};
		}
	
	 static Stream<VariantContext> streamVcf(final File f) {
		final VCFFileReader r = new VCFFileReader(f,false);
		final CloseableIterator<VariantContext> iter = r.iterator();
		return StreamSupport.stream(new IterableAdapter<VariantContext>(iter).spliterator(), false).onClose(()->{iter.close();r.close();});
		}
	 static Stream<VariantContext> streamJeterVcf() {
		return streamVcf(JETER_VCF);
		}
	
	 static VCFHeader getVcfHeader(final File f) {
		final VCFFileReader r = new VCFFileReader(f,false);
		final VCFHeader h = r.getFileHeader();
		r.close();
		return h;
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
    
    @Test
    public void testVcfRename() throws IOException {
    	File output = new File(TEST_RESULTS_DIR,"jeter.vcf");
    	File mapFile = new File(TEST_RESULTS_DIR,"jeter.tsv");
    	PrintWriter pw = new PrintWriter(mapFile);
    	
    	
    	BufferedReader br=new BufferedReader(new FileReader(TOY_FA+".fai"));
    	br.lines().map(L->L.split("\t")[0]).
    		map(L->L+"\tnew"+L).
    		forEach(L->pw.println(L));
    	br.close();
    	pw.flush();pw.close();
    	
        Assert.assertEquals(0,new ConvertVcfChromosomes().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--mapping",mapFile.getPath(),
        		TOY_VCF_GZ
        	}));
        Assert.assertTrue( getVcfHeader(output).getContigLines().stream().allMatch(C->C.getID().startsWith("new")));
        Assert.assertTrue(streamVcf(output).allMatch(V->V.getContig().startsWith("new")));
        Assert.assertTrue( output.delete());
        Assert.assertTrue( mapFile.delete());
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
        
        Assert.assertEquals(
        		Long.parseLong(IOUtil.slurp(output).trim()),
        		streamVcf(new File(vcfPath)).count()
        		);
        Assert.assertTrue( output.delete());
    	}
    
    @Test
    public void testBioAlcidaeJdkFasta() throws IOException{
    	final File fasta = new File(TEST_RESULTS_DIR,"jeter.fa");
		PrintWriter pw = new PrintWriter(fasta);
		int n = 100;
		for(int i=0;i< n;++i)
			{
			pw.println(">G"+i+"\nGAATTC\nGAATT");
			}
		pw.flush();
		pw.close();
		
    	File output = new File(TEST_RESULTS_DIR,"jeter.txt");
    
        Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-e","print(stream().count());",
        		fasta.getPath()
        	}));
        Assert.assertEquals( n ,  Integer.parseInt(IOUtil.slurp(output).trim()));
        Assert.assertTrue( output.delete());
        Assert.assertTrue( fasta.delete());
    	}
    
    @Test
    public void testBioAlcidaeJdkFastq() throws IOException{
    	final File fastq = new File(TEST_RESULTS_DIR,"jeter.fq");
		PrintWriter pw = new PrintWriter(fastq);
		int n = 100;
		for(int i=0;i< n;++i)
			{
			pw.println("@G"+i+"\nGAATTC\n+\n######");
			}
		pw.flush();
		pw.close();
		
    	File output = new File(TEST_RESULTS_DIR,"jeter.txt");
    
        Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-e","print(stream().count());",
        		fastq.getPath()
        	}));
        Assert.assertEquals( n ,  Integer.parseInt(IOUtil.slurp(output).trim()));
        Assert.assertTrue( output.delete());
        Assert.assertTrue( fastq.delete());
    	}
    
    @Test
    public void testBioAlcidaeJdkSam() throws IOException{
    	final File output = new File(TEST_RESULTS_DIR,"jeter.txt");
    
        Assert.assertEquals(0,new BioAlcidaeJdk().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-e","print(stream().count());",
        		TOY_BAM
        	}));
        Assert.assertEquals( 12 ,  Integer.parseInt(IOUtil.slurp(output).trim()));
        Assert.assertTrue( output.delete());
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
    
    @Test
    public void testSamJdk() throws IOException{
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.sam");
    	Assert.assertEquals(0,new SamJdk().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-e","return record.getStart()%2==0;",
        		TOY_BAM
        	}));
        Assert.assertTrue( tmp.exists());
    	}
    @Test
    public void testBamToRaster() throws IOException{
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.png");
    	Assert.assertEquals(0,new Bam2Raster().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-R",TOY_FA,
        		"-r","ref:1-100",
        		TOY_BAM
        	}));
        Assert.assertTrue( tmp.exists());
    	}
    @Test
    public void testBamToSql() throws IOException{
    	File tmp = new File(TEST_RESULTS_DIR,"jeter.sql");
    	Assert.assertEquals(0,new BamToSql().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-R",TOY_FA,
        		"-r","ref:1-100",
        		TOY_BAM
        	}));
        Assert.assertTrue( tmp.exists());
    	}
    @Test
    public void testLowResBamToRaster() throws IOException{
	File tmp = new File(TEST_RESULTS_DIR,"jeter.png");
	Assert.assertEquals(0,new LowResBam2Raster().instanceMain(new String[]{
    		"-o",tmp.getPath(),
    		"-R",TOY_FA,
    		"-r","ref:1-100",
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
    	final File tmp=new File(TEST_RESULTS_DIR,"tmp.ped.vcf");
        Assert.assertEquals(0,new VcfInjectPedigree().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"--pedigree",PED01,
        		VCF01
        	}));
        Assert.assertTrue(tmp.exists());
        Assert.assertTrue(
        		getVcfHeader(tmp).
        		getOtherHeaderLines()
        		.stream().
        		anyMatch(H->H.getKey().equals("Sample"))
        		);
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
    	if(getVcfHeader(new File(vcfPath)).getNGenotypeSamples()==0) return;
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
        		"--variantsWinCount","50",
        		"--variantsWinShift","50",
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
    
    public void testBiostar59647() throws IOException{    
    	File tmp = new File(TEST_RESULTS_DIR, "tmp.xml");
    	Assert.assertEquals(0,new Biostar59647().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"-R",TOY_FA,
        		TOY_BAM
        		}));
    	}
    @Test
    public void testVcfFilterSo() throws IOException{   
    	File output = new File(TEST_RESULTS_DIR, "jeter.filrerso.vcf");
    	final AnnPredictionParser parser = new AnnPredictionParserFactory().createDefaultParser();
    	final SequenceOntologyTree tree = SequenceOntologyTree.getInstance();
    	String acn = "SO:0001583";
    	final SequenceOntologyTree.Term term =tree.getTermByAcn(acn);
    	final Set<SequenceOntologyTree.Term> terms = term.getAllDescendants();
    	Assert.assertNotNull(term);
    	Assert.assertTrue(terms.size()>1);
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-A",acn,
        		VCF01
        		}));
    	streamVcf(output).forEach(V->{
    		//System.err.println(V.getAttribute("ANN")+" vs "+ terms);
    		Assert.assertTrue(parser.getPredictions(V).stream().
    			flatMap(P->P.getSOTerms().stream()).
				anyMatch(T->terms.contains(T)));
    	});
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-A",acn,
        		"--invert",
        		VCF01
        		}));
    	streamVcf(output).forEach(V->{
    		Assert.assertFalse(parser.getPredictions(V).stream().
    			flatMap(P->P.getSOTerms().stream()).
				anyMatch(T->terms.contains(T)));
    	});

    	
    	
    	Assert.assertEquals(0,new VcfFilterSequenceOntology().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-A",acn,
        		"--rmatt",
        		VCF01
        		}));
    	Assert.assertTrue(streamVcf(output).findAny().isPresent());
    	
    	Assert.assertTrue(output.delete());
    	}
    @Test
    public void testVcfRebase() throws IOException{    
    	Assert.assertEquals(0,new VcfRebase().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-R",TOY_FA,
        		TOY_VCF_GZ}));
    	}
    @Test
    public void testCaseCtrlCanvas() throws IOException{    
    	File tmp = new File(TEST_RESULTS_DIR, "tmp.png");

    	Assert.assertEquals(0,new CaseControlCanvas.Main().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"--pedigree",PED01,
        		VCF01}));
    	}
    @Test
    public void testVcfStats() throws IOException{    
    	final File tmp = new File(TEST_RESULTS_DIR, "jeter.zip");

    	Assert.assertEquals(0,new VcfStats().instanceMain(new String[]{
        		"-o",tmp.getPath(),
        		"--pedigree",PED01,
        		VCF01}));
    	Assert.assertTrue( tmp.exists());
    	}
    @Test
    public void testVcfMoveFiltersToInfo() throws IOException{    

    	Assert.assertEquals(0,new VcfMoveFiltersToInfo().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		VCF01}));
    	Assert.assertTrue( JETER_VCF.exists());
    	}
    @Test(groups={"burden"},dependsOnMethods={"testVcfInjectPed"})
    public void testVcfFilterNotInPedigree() throws IOException{    
    	final File input =new File(TEST_RESULTS_DIR,"tmp.ped.vcf");
    	Assert.assertEquals(0,new VcfFilterNotInPedigree().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		input.getPath()
        		}));
    	Assert.assertTrue( JETER_VCF.exists());
    	}
    @Test
    public void testVcfBigWig() throws IOException{    

    	Assert.assertEquals(0,new VCFBigWig().instanceMain(new String[]{
        		"-o",JETER_VCF.getPath(),
        		"-T","Uniqueness35bp",
        		"-B","src/test/resources/Uniqueness35bp.bigWig",
        		VCF01}));
    	Assert.assertTrue( JETER_VCF.exists());
    	Assert.assertTrue(streamVcf(JETER_VCF).
    			anyMatch(V->V.hasAttribute("Uniqueness35bp")
    			));
    	}
    

    @Test(groups={"burden"},dependsOnMethods={"testVcfInjectPed"})
    public void testVcfMultiToOneAllele() throws IOException{   
    	final File input =new File(TEST_RESULTS_DIR,"tmp.ped.vcf");
    	final File output =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
    	Assert.assertEquals(0,new VcfMultiToOneAllele().instanceMain(new String[]{
    			"--samples",
        		"-o",output.getPath(),
        		input.getPath()
        		}));
    	Assert.assertTrue( output.exists());
    	Assert.assertFalse(streamVcf(output).anyMatch(V->V.getAlternateAlleles().size()>1));
    	}
    
    @Test(groups={"burden"},dependsOnMethods={"testVcfMultiToOneAllele"})
    public void testVcfBurdenFilterExac() throws IOException{    
    		final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
    	
	    	Assert.assertEquals(0,new VcfBurdenFilterExac().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		"--discardNotInExac",
	        		"--exac","src/test/resources/ExAC.r1.sites.vep.vcf.gz",
	        		input.getPath()
	        		}));
	    	Assert.assertTrue( JETER_VCF.exists());
	    	}
    @Test
    public void testVcfGnomad() throws IOException{    
    		final File manifest = new File(TEST_RESULTS_DIR,"gnomad.manifest");
    		PrintWriter pw = new PrintWriter(manifest);
    		pw.println("exome\t*\tsrc/test/resources/gnomad.exomes.r2.0.1.sites.vcf.gz");
    		pw.println("genome\t1\tsrc/test/resources/gnomad.genomes.r2.0.1.sites.1.vcf.gz");
    		pw.flush();
    		pw.close();
    		
	    	Assert.assertEquals(0,new VcfGnomad().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		"-m",manifest.getPath(),
	        		VCF01}));
	    	Assert.assertTrue( JETER_VCF.delete());
	    	Assert.assertTrue( manifest.delete());
	    	}
    @Test(dependsOnMethods={"testVcfInjectPed"})
    public void testVcfFisherH() throws IOException{
    		final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
	    	Assert.assertEquals(0,new VcfBurdenFisherH().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		input.getPath()
	        		}));
	    	Assert.assertTrue( JETER_VCF.exists());
	    	}
    @Test(groups={"burden"},dependsOnMethods={"testVcfMultiToOneAllele"})
    public void testVcfFisherV() throws IOException{    
    		final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
	    	Assert.assertEquals(0,new VcfBurdenFisherV().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		input.getPath()
	        		}));
	    	Assert.assertTrue( JETER_VCF.exists());
	    	}
    @Test(groups={"burden"},dependsOnMethods={"testVcfMultiToOneAllele"})
    public void testVcfBurdenMAF() throws IOException{    
    		final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
	    	Assert.assertEquals(0,new VcfBurdenMAF().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		input.getPath()
	        		}));
	    	Assert.assertTrue( JETER_VCF.delete());
	    	}
    @Test(groups={"burden"},dependsOnMethods={"testVcfMultiToOneAllele"})
    public void testVcfBurdenRscriptV() throws IOException{    
			final File input =new File(TEST_RESULTS_DIR,"tmp.multi2oneallele.vcf");
			final File output =new File(TEST_RESULTS_DIR,"tmp.R");
	    	Assert.assertEquals(0,new VcfBurdenRscriptV().instanceMain(new String[]{
	        		"-o",output.getPath(),
	        		input.getPath()
	        		}));
	    	Assert.assertTrue( output.exists());
	    	}
    @Test
    public void testVcfSortOnInfo() throws IOException{    
	    	Assert.assertEquals(0,new SortVcfOnInfo().instanceMain(new String[]{
	        		"-o",JETER_VCF.getPath(),
	        		"-T","DP",
	        		TOY_VCF_GZ
	        		}));
	    	Assert.assertTrue( JETER_VCF.delete());
	    	}
    @Test
    public void testVcf2Xml() throws IOException{    
		final File output =new File(TEST_RESULTS_DIR,"vcf.xml");
    	Assert.assertEquals(0,new Vcf2Xml().instanceMain(new String[]{
        		"-o",output.getPath(),
        		TOY_VCF_GZ
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testBamToWig() throws IOException{    
		final File output =new File(TEST_RESULTS_DIR,"jeter.wig");
    	Assert.assertEquals(0,new Bam2Wig().instanceMain(new String[]{
        		"-o",output.getPath(),
        		TOY_BAM
        		}));
    	Assert.assertTrue( output.delete());
    	}
    
    @Test
    public void testVcfCompareCallers() throws IOException{    
		final File output =new File(TEST_RESULTS_DIR,"jeter.zip");
    	Assert.assertEquals(0,new VcfCompareCallers().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-vcf1","VCF1",
        		"-vcf2","VCF2",
        		"--prefix","pfx",
        		TOY_VCF_GZ,TOY_VCF_GZ
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testNgsFilesSummary() throws IOException{   
		final File input =new File(TEST_RESULTS_DIR,"jeter.path.txt");
    	PrintWriter pw=new PrintWriter(input);
    	pw.println(TOY_VCF_GZ);
    	pw.println(TOY_BAM);
    	pw.flush();
    	pw.close();
		final File output =new File(TEST_RESULTS_DIR,"jeter.txt");
    	Assert.assertEquals(0,new NgsFilesSummary().instanceMain(new String[]{
        		"-o",output.getPath(),
        		input.getPath()
        		}));
    	Assert.assertTrue( output.delete());
    	Assert.assertTrue( input.delete());
    	}
    @Test
    public void testFindAllCoverageAtPosition() throws IOException{   
		final File input =new File(TEST_RESULTS_DIR,"jeter.path.txt");
    	PrintWriter pw=new PrintWriter(input);
    	pw.println(TOY_BAM);
    	pw.flush();
    	pw.close();
		final File output =new File(TEST_RESULTS_DIR,"jeter.txt");
    	Assert.assertEquals(0,new FindAllCoverageAtPosition().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-p","ref2:14",
        		input.getPath()
        		}));
    	Assert.assertTrue( output.delete());
    	Assert.assertTrue( input.delete());
    	}
    @Test
    public void testFindAVariation() throws IOException{   
		final File input =new File(TEST_RESULTS_DIR,"jeter.path.txt");
    	PrintWriter pw=new PrintWriter(input);
    	pw.println(TOY_VCF_GZ);
    	pw.flush();
    	pw.close();
		final File output =new File(TEST_RESULTS_DIR,"jeter.txt");
    	Assert.assertEquals(0,new FindAVariation().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-p","ref2:14",
        		input.getPath()
        		}));
    	Assert.assertTrue( output.delete());
    	Assert.assertTrue( input.delete());
    	}
    @Test
    public void testSAM4WebLogo() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.txt");
    	Assert.assertEquals(0,new SAM4WebLogo().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-r","ref2:10-20",
        		TOY_BAM
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testMiniCaller() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new MiniCaller().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-R",TOY_FA,
        		TOY_BAM
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testVCFTrios() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VCFTrios().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-p",PED01,
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testVcfMultiToOne() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfMultiToOne().instanceMain(new String[]{
        		"-o",output.getPath(),
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testSortVcfOnInfo() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new SortVcfOnInfo().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"-T","AC",
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testNoCallToHomRef() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfNoCallToHomRef().instanceMain(new String[]{
        		"-o",output.getPath(),
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	}
    
    @Test
    public void testVcfAmalgation() throws IOException{   
    	final File tmpXml = File.createTempFile("_tmp", ".xml",TEST_RESULTS_DIR);
    	PrintWriter pw = new PrintWriter(tmpXml);
    	pw.println(
    		"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<config>\n"+
        		"<vcfhead><count>10</count></vcfhead>"+
        		"<vcffilterjdk><filter>ODD</filter><expression>return variant.getStart()%2==0;</expression></vcffilterjdk>"+
        		"<vcftail><count>5</count></vcftail>"+
    		"</config>"
    		);
    	pw.flush();
    	pw.close();
    	
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfXmlAmalgamation().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--xml",tmpXml.getPath(),
        		VCF01
        		}));
    	Assert.assertTrue( output.delete());
    	tmpXml.delete();
    	}
    @Test
    public void testVcfBurdenFilterGenes() throws IOException{   
    	final File tmp = File.createTempFile("_tmp", ".txt",TEST_RESULTS_DIR);
    	PrintWriter pw = new PrintWriter(tmp);
    	pw.println("ISG15");
    	pw.flush();
    	pw.close();
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VcfBurdenFilterGenes().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--genes",tmp.getPath(),
        		"--filter","ABCD",
        		"src/test/resources/ExAC.r1.sites.vep.vcf.gz"
        		}));
    	Assert.assertTrue( streamVcf(output).count()>0L);
    	Assert.assertTrue( output.delete());
    	tmp.delete();
    	}
    @Test
    public void testVCFStripAnnotations() throws IOException{   
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new VCFStripAnnotations().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"--exclude","INFO/CSQ,ID,FILTER/*",
        		"src/test/resources/ExAC.r1.sites.vep.vcf.gz"
        		}));
    	Assert.assertTrue( streamVcf(output).count()>0L);
    	Assert.assertTrue( streamVcf(output).filter(V->V.hasAttribute("CSQ")).count()==0L);
    	Assert.assertTrue( streamVcf(output).filter(V->V.hasID()).count()==0L);
    	Assert.assertTrue( streamVcf(output).filter(V->V.isFiltered()).count()==0L);
    	Assert.assertTrue( output.delete());
    	}
    
    @Test
    public void testGff2kg() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.kg");
    	Assert.assertEquals(0,new Gff2KnownGene().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"src/test/resources/gencode.v19.annotation.gff3"
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testFixVcfMissingGenotypes() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.vcf");
    	Assert.assertEquals(0,new FixVcfMissingGenotypes().instanceMain(new String[]{
        		"-B",TOY_BAM,
        		"-o",output.getPath(),
        		TOY_VCF_GZ
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testFastqShuffle() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.fq.gz");
    	Assert.assertEquals(0,new FastqShuffle().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"./src/test/resources/SAMPLE1_GATGAATC_L002_R1_001.fastq.gz",
        		"./src/test/resources/SAMPLE1_GATGAATC_L002_R2_001.fastq.gz"
        		}));
    	Assert.assertTrue( output.delete());
    	}
    @Test
    public void testPad() throws IOException {
		final File output =new File(TEST_RESULTS_DIR,"jeter.fq.gz");
    	Assert.assertEquals(0,new PadEmptyFastq().instanceMain(new String[]{
        		"-o",output.getPath(),
        		"./src/test/resources/SAMPLE1_GATGAATC_L002_R1_001.fastq.gz"
        		}));
    	Assert.assertTrue( output.delete());
    	}
	}
