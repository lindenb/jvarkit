/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.calling;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.tools.misc.ConcatSam;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.IntervalParser;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceContig;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceGenome;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceGenomeFactory;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/*
BEGIN_DOC


## Cited-In

  * "Direct Head-to-Head Evaluation of Recombinant Adeno-associated Viral Vectors Manufactured in Human versus Insect Cells". Kondratov & al. Molecular Therapy. [https://doi.org/10.1016/j.ymthe.2017.08.003](https://doi.org/10.1016/j.ymthe.2017.08.003).

## Example

```bash
$  java -jar dist/minicaller.jar -R ref.fa  bam.list > out.vcf
```

```
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Depth ReforAlt|Strand : RF,RR,AF,AR">
##FORMAT=<ID=DPG,Number=G,Type=Integer,Description="Depth for each allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Variant is indel">
##contig=<ID=rotavirus,length=1074>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	4	.	T	A	.	.	DP=65	GT:DP:DP4:DPG	0/1:22:20,0,2,0:20,2	./.	./.	./.
rotavirus	5	.	T	A	.	.	DP=83	GT:DP:DP4:DPG	0:27:27,0,0,0:27	./.	0/1:25:20,0,5,0:20,5	./.
rotavirus	6	.	T	A	.	.	DP=97	GT:DP:DP4:DPG	0:30:30,0,0,0:30	0/1:21:20,0,1,0:20,1	0/1:33:31,0,2,0:31,2	./.
rotavirus	7	.	T	A	.	.	DP=112	GT:DP:DP4:DPG	0/1:38:36,0,2,0:36,2	0/1:23:21,0,2,0:21,2	0:37:37,0,0,0:37	./.
rotavirus	8	.	A	C	.	.	DP=122	GT:DP:DP4:DPG	0/1:41:38,0,3,0:38,3	0/1:26:25,0,1,0:25,1	0/1:40:38,0,2,0:38,2	./.
rotavirus	9	.	A	C	.	.	DP=139	GT:DP:DP4:DPG	0/1:46:44,0,2,0:44,2	0/1:29:27,0,2,0:27,2	0/1:48:44,0,4,0:44,4	./.
```

END_DOC
 */
@Program(name="minicaller",
	description="Simple and Stupid Variant Caller designed for @AdrienLeger2",
	keywords={"bam","sam","calling","vcf"}
	)
public class MiniCaller extends Launcher
    {
	private static final Logger LOG = Logger.build(MiniCaller.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--mindepth"},description="Min depth")
	private int min_depth = 20 ;
    @Parameter(names={"-R","--reference"},
    		description=ReferenceGenomeFactory.OPT_DESCRIPTION_FILE_ONLY,
    		required=true)
    private File fastaFile = null;
	@Parameter(names={"--groupby"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition samRecordPartition = SAMRecordPartition.sample;
	@Parameter(names={"-f","--filter"},description="[20171130](replaced with jexl expression). "+SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter readFilter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-r","--region"},description=IntervalParser.OPT_DESC)
	private String rgnStr  = null;

	
	private SAMSequenceDictionary dictionary=null;
    private VariantContextWriter variantContextWriter = null;
    private ReferenceGenome referenceGenome=null;
    //private final Map<String,Integer> sample2index=new TreeMap<>();
  
    private final List<MyVariantContext> buffer=new ArrayList<>();
    private double min_fraction_alt=1.0/1000.0;

    
    private static class AlleleData
	    {
    	@SuppressWarnings("unused")
		final Allele alt;
    	 final int count_strands[]={0,0};
    	//final byte qual;
    	AlleleData(final Allele alt) {
    		this.alt = alt;
    		}
    	void incr(boolean negativeStrand)
    		{
    		this.count_strands[negativeStrand?1:0]++;
    		}
    	int count() {
    		return count_strands[0]+count_strands[1];
    		}
	    }
    
    private class SampleData
        {
    	final MyVariantContext owner;
    	@SuppressWarnings("unused")
		final String sampleName;
    	final Map<Allele,AlleleData> alleleMap=new HashMap<>();
    	SampleData(final MyVariantContext owner,final String sampleName) {
    		this.owner = owner;
    		this.sampleName = sampleName;
    		}
    	
    	public AlleleData getAllele(Allele alt)
    		{
    		if(owner.ref.equals(alt, true /* ignore state */))
    			{
    			alt = owner.ref;
    			}
    		AlleleData ad = this.alleleMap.get(alt);
    		if(ad==null)
    			{
    			ad = new AlleleData(alt);
    			this.alleleMap.put(alt, ad);
    			}
    		return ad;
    		}
    	
        }
    
    /** contains available information for everyone at contig=tid,pos0, ref */
    private class MyVariantContext
        implements Comparable<MyVariantContext>,
        Locatable
        {
    	/** contig ig */
        int tid;
        /** position 0 */
        int pos0;
        /** REF allele */
        Allele ref;
        /** everybody */
        final Map<String,SampleData> sampleData = new HashMap<>();

        @Override
        public String getContig()
        	{
        	return dictionary.getSequence(this.tid).getSequenceName();
        	}
        @Override
        public int getStart() {
        	return pos0+1;
        	}
        @Override
        public int getEnd() {
        	return pos0 + ref.length();
        	}
        
        @Override
        public int compareTo(final MyVariantContext o2)
            {
            int i = this.tid - o2.tid;
            if( i !=0) return i;
            i =this.pos0 - o2.pos0;
            if( i !=0) return i;
            i = this.ref.compareTo(o2.ref);
            return i;
            }
        
        public SampleData getSample(final String sampleName)
        	{
        	SampleData sd = this.sampleData.get(sampleName);
        	if(sd==null) {
        		sd = new SampleData(this,sampleName);
        		this.sampleData.put(sampleName,sd);
        		}
        	return sd;
        	}
        
        
        @Override
        public String toString() {
        	return getContig()+":"+pos0+" "+this.sampleData.size();
        	}
        
        void print()
        	{
        	VariantContext ctx=make();
        	if(ctx==null) return;
        	variantContextWriter.add(ctx);
        	}

        VariantContext make()
            {
        	boolean indel=this.ref.getBaseString().length()!=1;
            VariantContextBuilder vcb=new
                    VariantContextBuilder();
            vcb.chr(this.getContig());
            vcb.start(this.getStart());
            
            final List<Genotype> genotypes=new ArrayList<>();
            final Set<Allele> alleles=new TreeSet<Allele>();
            int total_depth=0;

            for(final String sampleName : this.sampleData.keySet())
                {
            	final SampleData sd = this.sampleData.get(sampleName);
            	final Counter<Allele> count_alleles = new Counter<Allele>();
                int dp4[]=new int[]{0,0,0,0};
            	
            	for(final Allele allele: sd.alleleMap.keySet())
            		{
            		if(allele.isNonReference() && allele.getDisplayString().equals("N")) continue;
            		alleles.add(allele);
            		final AlleleData ad = sd.alleleMap.get(allele);
            		count_alleles.incr(allele,ad.count());
            		if(allele.isReference())
            			{
            			dp4[0] += ad.count_strands[0];
                    	dp4[1] += ad.count_strands[1];
            			}
            		else
            			{
            			dp4[2] += ad.count_strands[0];
                    	dp4[3] += ad.count_strands[1];
            			}
            		
            		}
            	
                
            		
                
                total_depth+= count_alleles.getTotal();
                if(count_alleles.getTotal()> MiniCaller.this.min_depth)
                    {
                	final ArrayList<Allele> sample_alleles=new ArrayList<>(count_alleles.getCountCategories());
                	final ArrayList<Integer> sample_depths=new ArrayList<>(count_alleles.getCountCategories());
                	for(final Allele a: count_alleles.keySetDecreasing())
                		{
                		//skip if fraction of variant too low
                		if((float)count_alleles.count(a)/(float)count_alleles.getTotal() < MiniCaller.this.min_fraction_alt)
                			{
                			continue;
                			}
                		if(a.getBaseString().length()!=1) indel=true;
                		
                		sample_alleles.add(a);
                		sample_depths.add((int)count_alleles.count(a));
                		}
                	if(!sample_alleles.isEmpty())
	                	{
	                	final GenotypeBuilder gb=new GenotypeBuilder(sampleName, sample_alleles);
	                	gb.DP((int)count_alleles.getTotal());
	                	gb.attribute("DPG", sample_depths);
	                	gb.attribute("DP4",Arrays.asList(dp4));
	                	final Genotype gt = gb.make();
	                    alleles.addAll(sample_alleles);
	                    genotypes.add(gt);
	                	}
                	}
               
                }
            if(genotypes.isEmpty()) return null;
            
            alleles.add(this.ref);
            if(indel) vcb.attribute("INDEL", Boolean.TRUE);
            vcb.attribute("DP", total_depth);
            vcb.genotypes(genotypes);
            
            final List<Allele> orderedAlleles = new ArrayList<>(alleles);
            vcb.alleles(orderedAlleles);
            final int an = genotypes.stream().mapToInt(G->G.getAlleles().size()).sum();
            
            vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY, an);
            List<Integer> ac = new ArrayList<>();
            for(final Allele alt: orderedAlleles)
            	{
            	if(alt.isReference()) continue;
            	ac.add( (int) genotypes.stream().flatMap(G->G.getAlleles().stream())
            			.filter(A->A.equals(alt)).
            			count());
            	}
           
            if(!ac.isEmpty())
            	{
            	vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, ac);
            	vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,
            			ac.stream().
            				map(AC->AC/(double)an).
            				collect(Collectors.toList())
            				);
            	}
            
            
            /*int max_length=0;
            for(Allele a:alleles)
            	{
            	max_length=Math.max(a.getBaseString().length(), max_length);
            	}*/
            vcb.stop(this.pos0+this.ref.getBaseString().length());

            final VariantContext ctx= vcb.make();
            if(ctx.getAlternateAlleles().isEmpty()) return null;   
            return ctx;
            }

        }



    public MiniCaller() {
        }

    private MyVariantContext findContext(final int tid,final int pos0,final Allele ref)
        {
        int idx=0;
        for(idx=0; idx<this.buffer.size(); ++idx)
            {
            MyVariantContext ctx = this.buffer.get(idx);
            if(ctx.tid < tid) continue;
            if(ctx.tid > tid) break;
            if(ctx.pos0 < pos0) continue;
            if(ctx.pos0 > pos0) break;
            int i = ctx.ref.compareTo(ref);
            if(i< 0) continue;
            if(i>0 ) break;
            //System.err.println("ok got "+ctx);
            return ctx;
            }
        //System.err.println("new at "+tid+":"+pos0+" idx="+idx+" N="+this.buffer);
        final MyVariantContext ctx=new MyVariantContext();
        ctx.tid=tid;
        ctx.pos0=pos0;
        ctx.ref=ref;
        this.buffer.add(idx, ctx);
        return ctx;
        }
    
    @Override
    public int doWork(final List<String> args) {
    	
    	ConcatSam.ConcatSamIterator iter=null;
        try {
            
            if(this.fastaFile==null)
                {
            	LOG.error("no REF");
                return -1;
                }
            
          
            
            /* load faid */
           

            final ReferenceGenomeFactory referenceGenomeFactory = new ReferenceGenomeFactory();
            this.referenceGenome= referenceGenomeFactory.openFastaFile(this.fastaFile);
            this.dictionary = this.referenceGenome.getDictionary();
            if(this.dictionary==null) {
            	LOG.error(JvarkitException.FastaDictionaryMissing.getMessage(this.fastaFile.getPath()));
            	}
            
            /* create sam record iterator */
            
            iter = new  ConcatSam.Factory().
	            addInterval(this.rgnStr).
	            setEnableUnrollList(true).
	            open(args);
            
            final SAMFileHeader samFileheader= iter.getFileHeader();
            final SAMSequenceDictionary dict= samFileheader.getSequenceDictionary();
            if(dict==null)
            	{
            	LOG.error(JvarkitException.BamDictionaryMissing.getMessage(String.join(", ", args)));
            	return -1;
            	}
           
            
            if(!SequenceUtil.areSequenceDictionariesEqual(dict,this.dictionary))
            	{
            	LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict,this.dictionary));
            	return -1;
            	}
            
            final List<SAMReadGroupRecord> groups = samFileheader.getReadGroups();
            if(groups==null || groups.isEmpty())
            	{
            	LOG.error("No group defined in input");
            	return -1;
            	}
            final Set<String> sampleSet=groups.stream().
            		map(srgr->this.samRecordPartition.apply(srgr,samRecordPartition.name())).
            		collect(Collectors.toSet());
          
                


            /* create VCF metadata */
            final Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY));
            metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY));
            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));
            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY));
            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY));
            metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY));

            metaData.add(new VCFFormatHeaderLine(
                    "DPG",
                    VCFHeaderLineCount.G,//one value of each genotype
                    VCFHeaderLineType.Integer,
                    "Depth for each allele"));
           
           
            metaData.add(new VCFFormatHeaderLine(
                    "DP4",
                    4,
                    VCFHeaderLineType.Integer,
                    "Depth ReforAlt|Strand : RF,RR,AF,AR"));
            
           
            metaData.add(new VCFInfoHeaderLine(
                    "INDEL",
                    1,
                    VCFHeaderLineType.Flag,
                    "Variant is indel"));
            

   
            //addMetaData(metaData);
            
            final VCFHeader vcfHeader=new VCFHeader(
                    metaData , sampleSet
                    );
            vcfHeader.setSequenceDictionary(this.dictionary);
            /* create variant context */
            this.variantContextWriter = super.openVariantContextWriter(outputFile);
            this.variantContextWriter.writeHeader(vcfHeader);

            ReferenceContig genomicSeq=null;
            SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(this.dictionary);
            for(;;)
                {
                SAMRecord rec=null;
                
                
                if(iter.hasNext())
                    {
                    rec=progress.watch(iter.next());
                    if(rec.getReadUnmappedFlag()) continue;
                    if(this.readFilter.filterOut(rec)) continue;
                  
                    /* flush buffer if needed */
                    while(!this.buffer.isEmpty() &&
                    	(this.buffer.get(0).tid < rec.getReferenceIndex() ||
                    	(this.buffer.get(0).tid == rec.getReferenceIndex() && (this.buffer.get(0).getEnd()) < rec.getAlignmentStart())))
                    	{
                    	this.buffer.remove(0).print();
                    	}
                    /* get genomic sequence at this position */
                    if(genomicSeq==null ||
                            !genomicSeq.getContig().equals(rec.getContig()))
                            {
                            genomicSeq = this.referenceGenome.getContig(rec.getContig());
                            }
                    final Cigar cigar= rec.getCigar();
                    if(cigar==null) continue;
                    int readPos=0;
                    int refPos0 = rec.getAlignmentStart() -1;//0 based-reference
                    final byte bases[]=rec.getReadBases();
                    final byte quals[]=rec.getBaseQualities();
                    final String sampleName= this.samRecordPartition.getPartion(rec,samRecordPartition.name());
                    
                    for(final CigarElement ce: cigar.getCigarElements())
                        {
                        final CigarOperator op =ce.getOperator();
                        switch(op)
                            {
                            case P: break;
                            case H: break;
                            case S: readPos+=ce.getLength(); break;
                            case N://go
                            case D:
                                {
                                if(refPos0>0)// we need base before deletion
	                                {
                                	char refBase=genomicSeq.charAt(refPos0-1);/* we use base before deletion */
                                	final StringBuilder sb=new StringBuilder(ce.getLength());
	                            	sb.append(refBase);
	                                for(int i=0;i< ce.getLength();++i)
	                                	{
	                                	sb.append(genomicSeq.charAt(refPos0+i));
	                                	}
	                                findContext(
                                            rec.getReferenceIndex(),
                                            refPos0-1,//we use base *before deletion */
                                            Allele.create(sb.toString(), true)
                                            ).
		                                getSample(sampleName).
		                                getAllele(Allele.create(String.valueOf(refBase),false)).
		                                incr( rec.getReadNegativeStrandFlag());
	                                
	                                }
                                refPos0+= ce.getLength();
                                break;
                                }
                            case I:
                                {
                                if(refPos0>0)
	                                {
                                	//float qual=0;
                                	char refBase=Character.toUpperCase( genomicSeq.charAt(refPos0-1));
                                	final StringBuilder sb=new StringBuilder(1+ce.getLength());
	                                sb.append(refBase);
	                                for(int i=0;i< ce.getLength();++i)
	                                	{
	                                	sb.append((char)bases[readPos+i]);
	                                	//qual+=(readPos + i < quals.length?quals[ readPos + i]:0);
	                                	}
	                                findContext(
                                            rec.getReferenceIndex(),
                                            refPos0-1,//we use base *before deletion */
                                            Allele.create(String.valueOf(refBase), true)
                                            ).
	                                		getSample(sampleName).
	                                		getAllele(Allele.create(sb.toString().toUpperCase(),false)).
	                                		incr( rec.getReadNegativeStrandFlag())
	                                		;
	                                }
                                readPos+=ce.getLength();
                                break;
                                }
                            case EQ: case M: case X:
                                {
                                for(int i=0; i< ce.getLength();++i)
                                    {
                                    findContext(
                                            rec.getReferenceIndex(),
                                            refPos0 + i,
                                            Allele.create(String.valueOf(genomicSeq.charAt( refPos0 + i)), true)
                                    		).
                                		getSample(sampleName).
                                		getAllele(Allele.create(String.valueOf((char)bases[ readPos + i ]),false)).
                                    	incr( rec.getReadNegativeStrandFlag())
                                    	;
                                    }
                                readPos+=ce.getLength();
                                refPos0+= ce.getLength();
                                break;
                                }

                            default : throw new IllegalStateException("Case statement didn't deal with cigar op: "+ op);
                            }
                        }
                    }
                else
                    {
                    break;
                    }
                }
            
            while(!buffer.isEmpty()) buffer.remove(0).print();
            progress.finish();
            iter.close();iter=null;
            this.variantContextWriter.close();this.variantContextWriter=null;
            return RETURN_OK;
            }
        catch (Exception e)
            {
        	LOG.error(e);
            return -1;
            }
        finally
            {
        	 CloserUtil.close(iter);
            CloserUtil.close(this.referenceGenome);
            CloserUtil.close(this.variantContextWriter);
            }
        }

    /**
    * @param args
    */
    public static void main(String[] args)
        {
        new MiniCaller().instanceMainWithExit(args);
        }

    }
