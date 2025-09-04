package com.github.lindenb.jvarkit.tools.vcfisec;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;



public class VcfISecTest  {
	


private List<String> getIndexedesFiles() {
	final  TestSupport support =new TestSupport();
	final List<String> L=new ArrayList<>();
	L.add("test_vcf01.vcf");
	L.add("manta.B00GWGD.vcf.gz");
	L.add("rotavirus_rf.vcf.gz");
	L.add("S5.vcf.gz");
	L.add("S4.vcf.gz");
	L.add("toy.vcf.gz");
	L.add("toy.bcf");
	L.add("rotavirus_rf.ann.vcf.gz");
	L.add("manta.D000Q1R.vcf.gz");
	L.add("S2.vcf.gz");
	L.add("rotavirus_rf.freebayes.vcf.gz");
	L.add("ExAC.r1.sites.vep.vcf.gz");
	L.add("gnomad.genomes.r2.0.1.sites.1.vcf.gz");
	L.add("S3.vcf.gz");
	L.add("S1.vcf.gz");
	L.add("manta.B00GWIU.vcf.gz");
	L.add("gnomad.exomes.r2.0.1.sites.vcf.gz");
	L.add("gnomad_v2_sv.sites.vcf.gz");
	L.add("manta.B00I9CJ.vcf.gz");
	L.add("rotavirus_rf.unifiedgenotyper.vcf.gz");
	L.add("rotavirus_rf.unifiedgenotyper.vcf.gz");
	return L.stream().map(S->support.resource(S)).collect(Collectors.toList());
	}

@DataProvider(name = "all-indexed-vcf-files")
public Object[][] getData() {
	return getIndexedesFiles().
			stream().
			map(S->new Object[] {S}).toArray(N->new Object[N][]);
}

@Test(dataProvider="all-indexed-vcf-files")
public void doBasicTest(final String vcfpath) throws Exception
	{
	
	final  TestSupport support =new TestSupport();
	try 
		{
		Path vcfList = support.createTmpPath(".list");
		Files.writeString(vcfList, getIndexedesFiles().stream().filter(S->!S.endsWith(vcfpath)).collect(Collectors.joining("\n")));		
		Path outVcf = support.createTmpPath(".vcf");

		
		final List<String> args= new ArrayList<>();
		args.add("-o");
		args.add(outVcf.toString());
		args.add("--file-list");
		args.add(vcfList.toString());
		args.add(vcfpath);
		
		Assert.assertEquals(new VcfISec().instanceMain(args), 0);
		}
	finally
		{
		support.removeTmpFiles();
		}
	}


}
