package com.github.lindenb.jvarkit.tools.vcfphasedpep;

import java.io.IOException;
import org.testng.annotations.Test;
import java.nio.file.Path;
import org.testng.Assert;

import com.github.lindenb.jvarkit.tools.tests.TestSupport;

public class VcfPhasedPeptideTest {
    @Test
    public void basicTest() throws IOException
            {
    		final TestSupport support =new TestSupport();
            try {
                    final Path out = support.createTmpPath(".json");
                    Assert.assertEquals(new VcfPhasedPeptide().instanceMain(new String[] {
                            "-R",support.resource("rotavirus_rf.fa"),
                            "-o",out.toString(),
                            "--gtf",support.resource("rotavirus_rf.gtf.gz"),
                            support.resource("S1.vcf.gz")
                            }),0);
                    }
            finally
                    {
                    support.removeTmpFiles();
                    }
            }

}