package com.github.lindenb.jvarkit.util;

import com.github.lindenb.jvarkit.tools.misc.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.misc.VcfTail;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VCFFilterJS;
import com.github.lindenb.jvarkit.tools.vcfgo.VcfGeneOntology;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

/**
 * because I'm lazy: just a central point of compilation when
 * jvarkit used as a library
 * @author lindenb
 *
 */
class Library{
	
private Library()
	{
	@SuppressWarnings("unused")
	Class<?> clazz[]={
			VCFUtils.class,
			VepPredictionParser.class,
			SnpEffPredictionParser.class,
			VCFFilterJS.class,
			VcfGeneOntology.class,
			VcfFilterSequenceOntology.class,
			VcfHead.class,
			VcfTail.class
			};
	}
}
