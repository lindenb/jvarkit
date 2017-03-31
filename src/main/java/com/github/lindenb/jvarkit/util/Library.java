package com.github.lindenb.jvarkit.util;

import com.github.lindenb.jvarkit.tools.bioalcidae.BioAlcidae;
import com.github.lindenb.jvarkit.tools.groupbygene.GroupByGene;
import com.github.lindenb.jvarkit.tools.misc.BedIndexTabix;
import com.github.lindenb.jvarkit.tools.misc.VCFPolyX;
import com.github.lindenb.jvarkit.tools.misc.VCFShuffle;
import com.github.lindenb.jvarkit.tools.misc.VcfCadd;
import com.github.lindenb.jvarkit.tools.misc.VcfCutSamples;
import com.github.lindenb.jvarkit.tools.misc.VcfHead;
import com.github.lindenb.jvarkit.tools.misc.VcfIndexTabix;
import com.github.lindenb.jvarkit.tools.misc.VcfMultiToOneAllele;
import com.github.lindenb.jvarkit.tools.misc.VcfTail;
import com.github.lindenb.jvarkit.tools.onesamplevcf.VcfMultiToOne;
import com.github.lindenb.jvarkit.tools.sortvcfonref.SortVcfOnInfo;
import com.github.lindenb.jvarkit.tools.vcfbed.VCFBed;
import com.github.lindenb.jvarkit.tools.vcfbed.VCFBedSetFilter;
import com.github.lindenb.jvarkit.tools.vcfbigwig.VCFBigWig;
import com.github.lindenb.jvarkit.tools.vcfcmp.VcfCompareCallersOneSample;
import com.github.lindenb.jvarkit.tools.vcfcmp.VcfIn;
import com.github.lindenb.jvarkit.tools.vcfconcat.VcfConcat;
import com.github.lindenb.jvarkit.tools.vcffilterjs.VCFFilterJS;
import com.github.lindenb.jvarkit.tools.vcffilterso.VcfFilterSequenceOntology;
import com.github.lindenb.jvarkit.tools.vcfgo.VcfGeneOntology;
import com.github.lindenb.jvarkit.tools.vcfmerge.VCFMerge2;
import com.github.lindenb.jvarkit.tools.vcfstripannot.VCFReplaceTag;
import com.github.lindenb.jvarkit.tools.vcftrios.VCFTrios;
import com.github.lindenb.jvarkit.util.igv.IgvSocket;
import com.github.lindenb.jvarkit.util.vcf.IndexedVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.swing.DefaultVcfTable;
import com.github.lindenb.jvarkit.util.vcf.swing.InfoTreeModel;

/**
 * because I'm lazy: just a central point of compilation when
 * jvarkit used as a library
 * @author lindenb
 *
 */
class Library{
	
@Deprecated
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
			VcfTail.class,
			VCFMerge2.class,
			VcfIn.class,
			GroupByGene.class,
			VCFTrios.class,
			VcfCutSamples.class,
			VcfCompareCallersOneSample.class,
			VcfConcat.class,
			InfoTreeModel.class,
			DefaultVcfTable.class,
			VCFBigWig.class,
			VCFBed.class,
			VCFPolyX.class,
			VcfCadd.class,
			VCFShuffle.class,
			BioAlcidae.class,
			SortVcfOnInfo.class,
			VCFBedSetFilter.class,
			VCFReplaceTag.class,
			VcfIndexTabix.class,IndexedVcfFileReader.class,
			IgvSocket.class,
			VcfMultiToOneAllele.class,
			BedIndexTabix.class,
			VcfMultiToOne.class
			};
	}
}
