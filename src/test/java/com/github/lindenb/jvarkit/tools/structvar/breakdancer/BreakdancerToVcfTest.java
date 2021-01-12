package com.github.lindenb.jvarkit.tools.structvar.breakdancer;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;

import org.testng.Assert;
import org.testng.annotations.Test;

import com.github.lindenb.jvarkit.tests.AlsoTest;
import com.github.lindenb.jvarkit.tools.tests.TestSupport;
import com.github.lindenb.jvarkit.util.jcommander.LauncherTest;

@AlsoTest(LauncherTest.class)
public class BreakdancerToVcfTest {
	private final TestSupport support = new TestSupport();
	@Test
	public void test01() throws IOException
		{
		try {
			final Path dataFile=support.createTmpPath(".out");
			
			final PrintWriter pw = new PrintWriter(Files.newBufferedWriter(dataFile));
			pw.println("#Software: 1.4.1-unstable-10-fdfe9f2-dirty (commit fdfe9f2-dirty)");
			pw.println("#Command: bdfast -o21 inv_del_bam_config ");
			pw.println("#Library Statistics:");
			pw.println("#NA19238_chr21_del_inv.bam\tmean:475.76\tstd:28.67\tuppercutoff:532.53\tlowercutoff:311.36\treadlen:90\tlibrary:H_IJ-NA19238-NA19238-extlibs\treflen:5626088\tseqcov:0.0410481\tphycov:0.108495\t2:42\t3:10\t32:10");
			pw.println("#NA19240_chr21_del_inv.bam\tmean:467.59\tstd:31.91\tuppercutoff:525.5\tlowercutoff:277.03\treadlen:90\tlibrary:H_IJ-NA19240-NA19240-extlibs\treflen:5626088\tseqcov:0.0374168\tphycov:0.0971984\t1:2\t3:8\t8:2\t32:8");
			pw.println("#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib\tNA19238_chr21_del_inv.bam\tNA19240_chr21_del_inv.bam");
			pw.println("21\t29185056\t23+2-\t21\t29185377\t23+2-\tINS\t-226\t99\t2\tNA19238_chr21_del_inv.bam|1:NA19240_chr21_del_inv.bam|1\tNA\tNA");
			pw.println("21\t29185462\t23+2-\t21\t29186122\t2+22-\tDEL\t545\t99\t21\tNA19238_chr21_del_inv.bam|21\t176.58\t167.89");
			pw.println("21\t34807694\t8+6-\t21\t34808852\t8+6-\tINS\t-304\t99\t3\tNA19238_chr21_del_inv.bam|1:NA19240_chr21_del_inv.bam|2\tNA\tNA");
			pw.println("21\t34808937\t8+6-\t21\t34809799\t3+4-\tINV\t737\t99\t2\tNA19240_chr21_del_inv.bam|2\t847.39\t878.83");

			pw.flush();
			pw.close();
			
			
			final Path out = support.createTmpPath(".vcf");
			Assert.assertEquals(new BreakdancerToVcf().instanceMain(new String[] {
				"-o",out.toString(),
				"-R",support.resource("human_b37.dict"),
				dataFile.toString()}),0);
			support.assertIsVcf(out);
			}
		catch(final Throwable err) {
			support.removeTmpFiles();
			}
		}

}
