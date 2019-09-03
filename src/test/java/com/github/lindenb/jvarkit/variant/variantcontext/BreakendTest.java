package com.github.lindenb.jvarkit.variant.variantcontext;

import java.util.Optional;

import org.testng.Assert;
import org.testng.annotations.Test;

public class BreakendTest {
@Test
public void test01() {
	Optional<Breakend> ob=Breakend.parse("A");
	Assert.assertFalse(ob.isPresent());
	ob=Breakend.parse("<*>");
	Assert.assertFalse(ob.isPresent());
	ob=Breakend.parse("G[17:198982]");
	Assert.assertFalse(ob.isPresent());

	
	for(final String s: new String[]{"G]17:198982]","G[17:198982["}) {
		ob=Breakend.parse(s);
		Assert.assertTrue(ob.isPresent());
		Breakend b = ob.get();
		
		Assert.assertEquals(b.getContig(), "17");
		Assert.assertEquals(b.getStart(),198982);
		Assert.assertEquals(b.getEnd(),198982);
		Assert.assertEquals(b.getLeftSequence(),"G");
		Assert.assertEquals(b.getRightSequence(),"");
		}
	
	}
}
