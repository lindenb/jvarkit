package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.engine.filters.MalformedReadFilter;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.engine.walkers.Allows;
import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.Downsample;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.downsampling.DownsampleType;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Requires(value={})
@Allows(value={DataSource.READS, DataSource.REFERENCE})
@Reference(window=@Window(start=-50,stop=50))
@Downsample(by= DownsampleType.BY_SAMPLE, toCoverage=1000000)
@By(DataSource.REFERENCE)
public class SoftClipAnnotator 
	extends RodWalker<Integer, Integer>  implements TreeReducible<Integer> {

   @ArgumentCollection
   protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

   @Output(doc="File to which variants should be written")
   protected VariantContextWriter vcfWriter = null;
   
   @Input(doc="sam Readers. We don't use the standard GATK '-I' option because:  https://github.com/broadinstitute/gatk-protected/issues/891 ",fullName="bams",shortName="bams",required=false)
   protected List<File> samFilenames=new ArrayList<>();
   
   @Argument(doc="Extend read location. see https://github.com/broadinstitute/gatk-protected/issues/891",shortName="xclip",fullName="extendclip",required=false)
   protected int extend4clip=200;
   
   private List<SamReader> samReaders=new ArrayList<>();
   private Map<String,List<SamReader>> sample2bam=new HashMap<>();
   
   private final VCFFormatHeaderLine numClipInGenotypeFormatHeaderLine = new VCFFormatHeaderLine("NSOFTCLIP",1,VCFHeaderLineType.Integer,"Number of reads containing a soft clipped base at this positions");
   private final VCFFilterHeaderLine filterHaveClipInVariant = new VCFFilterHeaderLine("SOFT_CLIP_IN_VARIANT","All called and non-homref genotypes contain soft-clipped reads");
   
   public void initialize() {
	   final Set<VCFHeaderLine> hInfo = new HashSet<>();
       final List<String> rodName = Arrays.asList(variantCollection.variants.getName());
       final Set<String> samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), rodName);
	   for(final String sample:samples)
	   	{
		   this.sample2bam.put(sample, new ArrayList<>());   
	   	}
       
       for ( final VCFHeaderLine line : GATKVCFUtils.getHeaderFields(getToolkit(),rodName) ) {
           if ( isUniqueHeaderLine(line, hInfo) )
               hInfo.add(line);
	   		}
       VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
       VCFHeader header2=new VCFHeader(vcfHeader);
       header2.addMetaDataLine(this.numClipInGenotypeFormatHeaderLine);
       header2.addMetaDataLine(this.filterHaveClipInVariant);
       vcfWriter.writeHeader(header2);
       final SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
       for(final File samFilename:IOUtil.unrollFiles(this.samFilenames,".bam"))
       	{
    	logger.info("opening "+samFilename);
    	final SamReader r=srf.open(samFilename);
    	
    	final Set<String> sampleset= new HashSet<>();
    	for(final SAMReadGroupRecord g:r.getFileHeader().getReadGroups())
    		{
    		if(g.getSample()==null || !this.sample2bam.containsKey(g.getSample())) continue;
    		sampleset.add(g.getSample());
    		}
    	if(sampleset.isEmpty())
    		{
    		logger.info("no VCF sample in "+samFilename);
    		CloserUtil.close(r);
    		continue;
    		}
    	this.samReaders.add(r);
    	
    	for(final String sample:sampleset)
    		{
    		if(!this.sample2bam.containsKey(sample)) continue;
    		this.sample2bam.get(sample).add(r);
    		}
       	}
       
       
   	}
   
   @Override
   public boolean includeReadsWithDeletionAtLoci() { return true; }
   
   public Integer map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {
	   if(tracker==null) return 0;
	   final Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
	   if ( VCs.isEmpty() )
	   	    {
		    return 0;
	   		}
	   
	   final Collection<ReadFilter> readFilters = super.getToolkit().getFilters();
	   
	  
       for ( final VariantContext ctx : VCs )
       	   {
    	   int observed_genotypes=0;
    	   int num_some_clipped_genotypes=0;
    	   final List<Genotype> genotypes=new ArrayList<>(ctx.getNSamples());
    	   for(int i=0;  i< ctx.getNSamples();++i)
    	   	{
    		final Genotype g=ctx.getGenotype(i);
    		if(!g.isCalled() || g.isNoCall() || g.isHomRef()) {
    			genotypes.add(g);
    			continue;
    		}
    		
    		
    		final List<SamReader> bamReaderForSample = this.sample2bam.get(g.getSampleName());
    		if(bamReaderForSample.isEmpty())
    			{
    			genotypes.add(g);
    			continue;
    			}
    		
    		observed_genotypes++;
    		
    		int numberOfClipsOverlapping=0;
    		
    		for(final SamReader samReader: bamReaderForSample)
    			{
    			final SAMRecordIterator iter= samReader.query(
    					ctx.getContig(),
    					Math.max(0, ctx.getStart()-extend4clip),
    					ctx.getEnd()+extend4clip
    					, false);
    			
    			
    			while(iter.hasNext())
    				{
    				
    				final SAMRecord samRecord = iter.next();
        			if(samRecord.getReadUnmappedFlag() || samRecord.getCigar()==null ||
        					!samRecord.getContig().equals(ctx.getContig()) ||
        					samRecord.getUnclippedEnd() < ctx.getStart() ||
        					samRecord.getUnclippedStart() > ctx.getEnd() ||
        					samRecord.getReadGroup()==null ||
        					!g.getSampleName().equals(samRecord.getReadGroup().getSample())
        					) continue;
        			
        			boolean ok=true;
        			
        			for(final ReadFilter readFilter:readFilters)
        				{
        				//this one cannot't be correctly initialized...
        				if(readFilter instanceof MalformedReadFilter) continue;
        				
        				if(readFilter.filterOut(samRecord))
        					{
        					ok=false;
        					break;
        					}
        				}
        			if(!ok) continue;
        			int refPos=samRecord.getUnclippedStart();
        			for(final CigarElement ce:samRecord.getCigar())
	    				{
						if( ce.getOperator().consumesReferenceBases() ||
							ce.getOperator().isClipping())
							{
							final int refEnd = refPos+ce.getLength();
							if(!(ctx.getEnd() < refPos || refEnd < ctx.getStart()))
								{
								//System.err.println("IN "+ce+" "+ctx.getStart()+"-"+ctx.getEnd()+" : "+refPos+"-"+refPos+ce.getLength());
		    					if(ce.getOperator().equals(CigarOperator.S))
		    						{
		    						++numberOfClipsOverlapping;
		    						}
		    					}
							refPos+=ce.getLength();
							}
	    				if(refPos> ctx.getEnd()) break;
	    				}    				
        			}
    			
    			iter.close();
    			}//end of loop over SamRecord
    		
    		if(numberOfClipsOverlapping==0)
    			{
    			genotypes.add(g);
    			}
    		else
    			{
    			GenotypeBuilder gb=new GenotypeBuilder(g);
    			gb.attribute(numClipInGenotypeFormatHeaderLine.getID(), numberOfClipsOverlapping);
    			genotypes.add(gb.make());
    			num_some_clipped_genotypes++;
    			}
    	   	}/* end loop oversam Reader */
    	 
    	 
    	 if(num_some_clipped_genotypes==0)
    	 	{
    		 vcfWriter.add(ctx);
    	 	}
    	 else
	    	 {
        	 final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
        	 vcb.genotypes(genotypes);
        	 if(observed_genotypes==num_some_clipped_genotypes)
	        	 {
	        	 vcb.filter(this.filterHaveClipInVariant.getID());
	        	 }
        	 vcfWriter.add(vcb.make());
	    	 }
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
	   for(final SamReader samReader:this.samReaders)
	   	{
		CloserUtil.close(samReader);   
	   	}
       logger.info("Processed " + result + " loci.\n");
}
   
   
   private static boolean isUniqueHeaderLine(VCFHeaderLine line, Set<VCFHeaderLine> currentSet) {
       if ( !(line instanceof VCFCompoundHeaderLine) )
           return true;

       for ( VCFHeaderLine hLine : currentSet ) {
           if ( hLine instanceof VCFCompoundHeaderLine && ((VCFCompoundHeaderLine)line).sameLineTypeAndName((VCFCompoundHeaderLine)hLine) )
               return false;
       }

       return true;
   	}
   
   
}
