package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.walkers.Allows;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Requires(value={})
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@Reference(window=@Window(start=-50,stop=50))
@Downsample(by= DownsampleType.BY_SAMPLE, toCoverage=1000000)
@By(DataSource.REFERENCE)
public class AltVariantAnnotator 
	extends RodWalker<Integer, Integer>  implements TreeReducible<Integer> {

   @ArgumentCollection
   protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

   @Output(doc="File to which variants should be written")
   protected VariantContextWriter vcfWriter = null;
   
   public void initialize() {
	   final Set<VCFHeaderLine> hInfo = new HashSet<>();
       final List<String> rodName = Arrays.asList(variantCollection.variants.getName());
       final Set<String> samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), rodName);
	   for ( final VCFHeaderLine line : GATKVCFUtils.getHeaderFields(getToolkit(),rodName) ) {
           if ( isUniqueHeaderLine(line, hInfo) )
               hInfo.add(line);
	   		}
       VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
       vcfWriter.writeHeader(vcfHeader);
   	}
   
   @Override
   public boolean includeReadsWithDeletionAtLoci() { return true; }
   
   public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
	   if(tracker==null) return 0;
	   Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
	   if ( VCs.isEmpty() )
	   	    {
		   return 0;
	   		}
	   
	   final ReadBackedPileup rbp = context.getBasePileup();
	   
	   
       for ( final VariantContext ctx : VCs )
       	{
    	  
    	   for(int i=0;rbp!=null && i< ctx.getNSamples();++i)
    	   	{
    		final Genotype g=ctx.getGenotype(i);
    		if(g.isNoCall() || g.isHomRef()) continue;
    		final ReadBackedPileup sampleReadBackedPileup =rbp.getPileupForSample(g.getSampleName());
    		 if(sampleReadBackedPileup==null || sampleReadBackedPileup.isEmpty()) 
    		 	{
    			//NO Read for this sample
    			continue;
    		 	}
    		 for(final GATKSAMRecord gatksamRecord:sampleReadBackedPileup.getReads())
    		 	{
    			boolean debug=(gatksamRecord.getReadName().equals("HWI-1KL149:89:HA4T6ADXX:1:1112:16630:78937") ||
    				gatksamRecord.getReadName().equals("HWI-1KL149:89:HA4T6ADXX:1:1112:16630:78937"));
    			if(debug)
    			{
    				System.err.println(gatksamRecord+" "+gatksamRecord.getCigarString());
    			}
    				
    			if(gatksamRecord.getReadUnmappedFlag() || gatksamRecord.getCigar()==null) continue;
    			int refPos=gatksamRecord.getUnclippedStart();
    			for(final CigarElement ce:gatksamRecord.getCigar())
    				{
    				if(debug) System.err.println(ce+" "+refPos+" "+gatksamRecord.getReadGroup().getSample());
    				if(refPos<= ctx.getStart() && ctx.getEnd()<=refPos+ce.getLength())
    					{
    					if(debug) System.err.println("IN "+ce+" "+ctx.getStart()+"-"+ctx.getEnd()+" : "+refPos+"-"+refPos+ce.getLength());
    					if(ce.getOperator().equals(CigarOperator.S))
    						{
    						System.err.println("YYYYes");
    						}
    					else if(ce.getOperator().isAlignment())
    						{
    						
    						}
    					}
    				
					if(ce.getOperator().consumesReferenceBases() ||ce.getOperator().isClipping())
						{
						refPos+=ce.getLength();
						}
    				if(refPos> ctx.getEnd()) break;
    				}
    			if(debug) System.exit(0);
    		 	}
    		
    	   	}
    	   
    	  vcfWriter.add(ctx);
       	}
	  
	   
	   return VCs.size();
   }
   @Override
   public Integer reduceInit() { return 0; }

   @Override
   public Integer reduce(Integer value, Integer sum) { return value + sum; }

   @Override
   public Integer treeReduce(Integer lhs, Integer rhs) {
       return lhs + rhs;
   }

   /**
    * Tell the user the number of loci processed and close out the new variants file.
    *
    * @param result  the number of loci seen.
    */
   public void onTraversalDone(Integer result) {
       logger.info("Processed " + result + " loci.\n");
}
   
   
   public static boolean isUniqueHeaderLine(VCFHeaderLine line, Set<VCFHeaderLine> currentSet) {
       if ( !(line instanceof VCFCompoundHeaderLine) )
           return true;

       for ( VCFHeaderLine hLine : currentSet ) {
           if ( hLine instanceof VCFCompoundHeaderLine && ((VCFCompoundHeaderLine)line).sameLineTypeAndName((VCFCompoundHeaderLine)hLine) )
               return false;
       }

       return true;
   	}
   
   private  interface AltAnnotation
  	{
	   public String getName();
	   public String getDescription();
	   public List<VCFInfoHeaderLine> getDescriptions();
  	}

   
   private abstract interface AltVariantAnnotation extends AltAnnotation
   	{
	   public VariantContext annotate(VariantContext ctx, ReferenceContext ref, AlignmentContext context); 
   	}
   private abstract interface AltGenotypeAnnotation extends AltAnnotation
  	{
	   public Genotype annotate(VariantContext ctx, Genotype genotype, ReferenceContext ref, AlignmentContext context); 
  	}
   
}
