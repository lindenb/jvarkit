package com.github.lindenb.jvarkit.bio;

import org.testng.Assert;
import org.testng.annotations.Test;

public class AcidNucleicsTest {
@Test
public void testAcidNucleic() {
	Assert.assertTrue( AcidNucleics.isATGC("ATCGATCGTA"));
	Assert.assertFalse(AcidNucleics.isATGC("ATCGATCGTAN"));
	Assert.assertTrue(AcidNucleics.isATGCN("ATCGATCGTAN"));
	Assert.assertFalse(AcidNucleics.isATGCN("ATCGATCGTANY"));
	Assert.assertTrue(AcidNucleics.isIUPAC("ATCGATCGTANY"));
	Assert.assertFalse(AcidNucleics.isIUPAC("ATCGATCGTANY#"));
	Assert.assertEquals(AcidNucleics.reverseComplement("AATTCG"),"CGAATT");
	Assert.assertEquals(AcidNucleics.degenerateToBases('N').length,4);
	Assert.assertEquals(AcidNucleics.degenerateToBases('W').length,2);
	}
}
