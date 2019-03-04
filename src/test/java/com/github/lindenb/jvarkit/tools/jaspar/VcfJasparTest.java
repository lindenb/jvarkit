package com.github.lindenb.jvarkit.tools.jaspar;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfJasparTest  {
	private final TestSupport support = new TestSupport();

	@DataProvider(name = "src1")
	public Object[][] createData1() {
		return new Object[][]{
			{support.resource("toy.vcf.gz"),support.resource("toy.fa")},
			{support.resource("rotavirus_rf.vcf.gz"),support.resource("rotavirus_rf.fa")}
			};
	}

	
	@Test(dataProvider="src1")
	public void test01(final String vcfin,final String ref) throws IOException {
		try {
    	final Path output = support.createTmpPath(".vcf");
    	final Path matrixFile = support.createTmpPath(".tmp");
    	final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(matrixFile));
    	pw.println(">MA0002.2\tRUNX1");
    	pw.println("A  [ 287  234  123   57    0   87    0   17   10  131  500 ]");
    	pw.println("C  [ 496  485 1072    0   75  127    0   42  400  463  158 ]");
    	pw.println("G  [ 696  467  149    7 1872   70 1987 1848  251   81  289 ]");
    	pw.println("T  [ 521  814  656 1936   53 1716   13   93 1339 1325 1053 ]");
    	pw.println("");
    	pw.println(">MA0003.3\tTFAP2A");
    	pw.println("A  [1706  137    0    0   33  575 3640 1012    0   31 1865 ]");
    	pw.println("C  [1939  968 5309 5309 1646 2682  995  224   31 4726  798 ]");
    	pw.println("G  [ 277 4340  139   11  658 1613  618 5309 5309  582 1295 ]");
    	pw.println("T  [1386   47    0  281 2972  438   56    0    0   21 1350 ]");
    	pw.println("");
    	pw.println(">MA0004.1\tArnt");
    	pw.println("A  [ 4 19  0  0  0  0 ]");
    	pw.println("C  [16  0 20  0  0  0 ]");
    	pw.println("G  [ 0  1  0 20  0 20 ]");
    	pw.println("T  [ 0  0  0  0 20  0 ]");
    	pw.println("");
    	pw.flush();
    	pw.close();
    	
        Assert.assertEquals(new VcfJaspar().instanceMain(new String[]{
        		"-o",output.toString(),
        		"-R",ref,
        		"-J",matrixFile.toString(),
        		vcfin
        	}),0);
        support.assertIsVcf(output);
	} finally {
		support.removeTmpFiles();
	}
    	}
}
