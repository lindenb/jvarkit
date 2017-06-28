/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.samfilter.SamFilterParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.semontology.Term;

/*
BEGIN_DOC


## Example

```bash
$  java -jar dist/minicaller.jar -R ref.fa  bam.list > out.vcf

```

END_DOC
 */
@Program(name="minicaller",
	description="Simple and Stupid Variant Caller designed for @AdrienLeger2",
	terms=Term.ID_0000005,
	keywords={"bam","sam","calling","vcf"}
	)
public class MiniCaller extends Launcher
    {
	private static final Logger LOG = Logger.build(MiniCaller.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--mindepth"},description="Min depth")
	private int min_depth = 20 ;
    @Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
    private File fastaFile = null;
	@Parameter(names={"--groupby"},description="Group Reads by")
	private SAMRecordPartition samRecordPartition = SAMRecordPartition.sample;
	@Parameter(names={"-f","--filter"},description=SamFilterParser.FILTER_DESCRIPTION,converter=SamFilterParser.StringConverter.class)
	private SamRecordFilter readFilter  = SamFilterParser.buildDefault();

	
	private SAMSequenceDictionary dictionary=null;
    private VariantContextWriter variantContextWriter = null;
    private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
    private Map<String,Integer> sample2index=new TreeMap<>();
    private List<String> samples=new ArrayList<>();
    private List<MyVariantContext> buffer=new ArrayList<>();
    private double min_fraction_alt=1.0/1000.0;

    
    private class Base2
        {
        Allele alt;
        @SuppressWarnings("unused")
		byte qual;
        int sample_index;
        int count_strands[]={0,0};
        Base2(Allele alt,byte qual,int sample_index,boolean negativeStrand)
            {
            this.alt=alt;
            this.qual=qual;
            this.sample_index=sample_index;
            this.count_strands[negativeStrand?1:0]=1;
            }
        public int count()
        	{
        	return this.count_strands[0]+this.count_strands[1];
        	}
        }
    private class MyVariantContext
        implements Comparable<MyVariantContext>
        {
        int tid;
        int pos0;
        Allele ref;
        List<Base2> bases = new ArrayList<Base2>();

        
        public String getReferenceName()
        	{
        	return dictionary.getSequence(this.tid).getSequenceName();
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
        
        void add(Base2 base)
            {
            if( base.alt.getBaseString().equals( this.ref.getBaseString() ))
                {
                base.alt = this.ref;
                }
            if(base.alt.getBaseString().equals("N")) return;
            
            //ugly, use std::lower_bound
            int idx=0;
            for(idx=0; idx<this.bases.size(); ++idx)
                {
                Base2 b2 = this.bases.get(idx);
                if(b2.sample_index < base.sample_index) continue;
                if(b2.sample_index > base.sample_index) break;
                int i = b2.alt.compareTo(base.alt);
                if(i< 0) continue;
                if(i>0 ) break;
                //go it
                b2.count_strands[0]+=base.count_strands[0];
                b2.count_strands[1]+=base.count_strands[1];
                return ;
                }
            //System.err.println("new at "+tid+":"+pos0+" idx="+idx+" N="+this.buffer);
            this.bases.add(idx, base);
            }
        
        @Override
        public String toString() {
        	return getReferenceName()+":"+pos0+" "+this.bases.size();
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
            vcb.chr(this.getReferenceName());
            vcb.start(this.pos0+1);
            
            List<Genotype> genotypes=new ArrayList<>();
            Set<Allele> alleles=new TreeSet<Allele>();
            int total_depth=0;

            for(int sample_index=0;sample_index <
                    samples.size();++sample_index)
                {
                Counter<Allele> count_alleles = new Counter<Allele>();
                int dp4[]=new int[]{0,0,0,0};

                for(Base2 v2: this.bases)
                    {
                    if(v2.sample_index!=sample_index) continue;
                    count_alleles.incr(v2.alt,v2.count());
                    
                    
                    
                    //cal dp4 RF,RR,AF,AR
                    if(v2.alt.isReference())
                    	{
                    	dp4[0] += v2.count_strands[0];
                    	dp4[1] += v2.count_strands[1];
                    	}
                    else
                    	{
                    	dp4[2] += v2.count_strands[0];
                    	dp4[3] += v2.count_strands[1];
                    	}
                    }
                
                total_depth+= count_alleles.getTotal();
                if(count_alleles.getTotal()> MiniCaller.this.min_depth)
                    {
                	ArrayList<Allele> sample_alleles=new ArrayList<>(count_alleles.getCountCategories());
                	ArrayList<Integer> sample_depths=new ArrayList<>(count_alleles.getCountCategories());
                	for(Allele a: count_alleles.keySetDecreasing())
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
	                	GenotypeBuilder gb=new GenotypeBuilder(
	                			MiniCaller.this.samples.get(sample_index),
	                			sample_alleles);
	                	gb.DP((int)count_alleles.getTotal());
	                	gb.attribute("DPG", sample_depths);
	                	gb.attribute("DP4",Arrays.asList(dp4));
	                	Genotype gt = gb.make();
	                    alleles.addAll(sample_alleles);
	                    genotypes.add(gt);
	                	}
                	}
               
                }
           
            
            alleles.add(this.ref);
            if(indel) vcb.attribute("INDEL", Boolean.TRUE);
            vcb.attribute("DP", total_depth);
            vcb.genotypes(genotypes);
            vcb.alleles(alleles);
            
            /*int max_length=0;
            for(Allele a:alleles)
            	{
            	max_length=Math.max(a.getBaseString().length(), max_length);
            	}*/
            vcb.stop(this.pos0+this.ref.getBaseString().length());

            VariantContext ctx= vcb.make();
            if(ctx.getAlternateAlleles().isEmpty()) return null;   
            return ctx;
            }

        }



    private MiniCaller()
        {

        }

    private MyVariantContext findContext(int tid,int pos0,Allele ref)
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
        MyVariantContext ctx=new MyVariantContext();
        ctx.tid=tid;
        ctx.pos0=pos0;
        ctx.ref=ref;
        this.buffer.add(idx, ctx);
        return ctx;
        }
    
    @Override
    public int doWork(final List<String> args) {
    	final Set<File> bamFileSet=new HashSet<File>();
        List<SamReader> readers = new ArrayList<SamReader>();
        try {
            for(String filename:args)
                {
                if(filename.endsWith(".list"))
                    {
                    BufferedReader
                    in=IOUtils.openFileForBufferedReading(new File(filename));
                    String line;
                    while((line=in.readLine())!=null)
                        {
                        if(line.trim().isEmpty()) continue;
                        bamFileSet.add(new File(line));
                        }
                    in.close();
                    }
                else
                    {
                    bamFileSet.add(new File(filename));
                    }
                }
            if(this.fastaFile==null)
                {
            	LOG.error("no REF");
                return -1;
                }
            if(bamFileSet.isEmpty())
	            {
		        LOG.error("No Bam Files");
		        return -1;
	            }
            /* load faid */
            this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(fastaFile);


            
            /** create merged SAMReader */
            List<File> bamFiles=new ArrayList<File>(bamFileSet);
            
            /* all readers */
            
            List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
            SamReaderFactory srf=SamReaderFactory.make();
            srf.validationStringency(ValidationStringency.LENIENT);
            
            /* open each bam file */
            for(File bamFile:bamFiles)
                {
                LOG.info("Opening "+bamFile);
                SamReader samReader = srf.open(bamFile);
                /* get header; check group and dict */
                SAMFileHeader header= samReader.getFileHeader();
                
                SAMSequenceDictionary dict= header.getSequenceDictionary();
                if(dict==null)
                	{
                	LOG.error("No SAMSequenceDictionary defined in "+bamFile);
                	return -1;
                	}
                if(this.dictionary==null)
                	{
                	this.dictionary=dict;
                	}
                
                if(!SequenceUtil.areSequenceDictionariesEqual(dict,this.dictionary))
                	{
                	LOG.error("Not same SAMSequenceDictionaries last was "+bamFile);
                	return -1;
                	}
                
                
                	
                List<SAMReadGroupRecord> groups = header.getReadGroups();
                if(groups==null || groups.isEmpty())
                	{
                	LOG.error("No group defined in "+bamFile);
                	return -1;
                	}
                
                for(final SAMReadGroupRecord srgr : groups)
                    {
                    String sampleName=this.samRecordPartition.apply(srgr);
                    if( sampleName == null ) sampleName=samRecordPartition.name();
                    if(!this.sample2index.containsKey(sampleName))
                        {
                        int sample_index=this.samples.size();
                        this.sample2index.put(sampleName, sample_index);
                        this.samples.add(sampleName);
                        }
                    }
                
                headers.add( header);
                readers.add(samReader);
                }

            /* create merged sam header */
            SamFileHeaderMerger merger=new SamFileHeaderMerger(
            		SAMFileHeader.SortOrder.coordinate,
            		headers,
            		false
            		);

            /* create sam record iterator */
            MergingSamRecordIterator iter= new MergingSamRecordIterator(
                    merger,
                    readers,
                    true
                    );

            /* create VCF metadata */
            Set<VCFHeaderLine> metaData=new HashSet<VCFHeaderLine>();
            metaData.add(new VCFFormatHeaderLine(
                    "GT",
                    1,
                    VCFHeaderLineType.String,
                    "Genotype"));
            metaData.add(new VCFFormatHeaderLine(
                    "DPG",
                    VCFHeaderLineCount.G,//one value of each genotype
                    VCFHeaderLineType.Integer,
                    "Depth for each allele"));
           
            metaData.add(new VCFFormatHeaderLine(
                    "DP",
                    1,
                    VCFHeaderLineType.Integer,
                    "Depth"));
           
            metaData.add(new VCFFormatHeaderLine(
                    "DP4",
                    4,
                    VCFHeaderLineType.Integer,
                    "Depth ReforAlt|Strand : RF,RR,AF,AR"));
            
            metaData.add(new VCFInfoHeaderLine(
                    "DP",
                    1,
                    VCFHeaderLineType.Integer,
                    "Depth"));
            
            metaData.add(new VCFInfoHeaderLine(
                    "INDEL",
                    1,
                    VCFHeaderLineType.Flag,
                    "Variant is indel"));
            

            /* add dict */
            metaData.addAll(VCFUtils.samSequenceDictToVCFContigHeaderLine(
            		this.dictionary
            		));
            //addMetaData(metaData);
            
            final VCFHeader vcfHeader=new VCFHeader(
                    metaData , this.sample2index.keySet()
                    );
            
            /* create variant context */
            this.variantContextWriter = super.openVariantContextWriter(outputFile);
            this.variantContextWriter.writeHeader(vcfHeader);

            GenomicSequence genomicSeq=null;
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
                    	(this.buffer.get(0).tid == rec.getReferenceIndex() && (this.buffer.get(0).pos0+1) < rec.getAlignmentStart())))
                    	{
                    	this.buffer.remove(0).print();
                    	}
                    /* get genomic sequence at this position */
                    if(genomicSeq==null ||
                            !genomicSeq.getChrom().equals(rec.getReferenceName()))
                            {
                            genomicSeq = new GenomicSequence(
                                    this.indexedFastaSequenceFile,
                                    rec.getReferenceName());
                            }
                    Cigar cigar= rec.getCigar();
                    int readPos=0;
                    int refPos0 = rec.getAlignmentStart() -1;//0 based-reference
                    byte bases[]=rec.getReadBases();
                    byte quals[]=rec.getBaseQualities();
                    String sampleName= this.samRecordPartition.getPartion(rec);
                    if( sampleName == null ) sampleName=samRecordPartition.name();
                    
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
	                            	StringBuilder sb=new StringBuilder(ce.getLength());
	                            	sb.append(refBase);
	                                for(int i=0;i< ce.getLength();++i)
	                                	{
	                                	sb.append(genomicSeq.charAt(refPos0+i));
	                                	}
	                                MyVariantContext v=findContext(
                                            rec.getReferenceIndex(),
                                            refPos0-1,//we use base *before deletion */
                                            Allele.create(sb.toString(), true)
                                            );
	                                v.add(new Base2(
                                            Allele.create(String.valueOf(refBase),false),
                                            (byte)0,
                                            this.sample2index.get(sampleName),
                                            rec.getReadNegativeStrandFlag()
                                            ));
	                                }
                                refPos0+= ce.getLength();
                                break;
                                }
                            case I:
                                {
                                if(refPos0>0)
	                                {
                                	float qual=0;
                                	char refBase=Character.toUpperCase( genomicSeq.charAt(refPos0-1));
	                                StringBuilder sb=new StringBuilder(1+ce.getLength());
	                                sb.append(refBase);
	                                for(int i=0;i< ce.getLength();++i)
	                                	{
	                                	sb.append((char)bases[readPos+i]);
	                                	qual+=quals[readPos+i];
	                                	}
	                                MyVariantContext v=findContext(
                                            rec.getReferenceIndex(),
                                            refPos0-1,//we use base *before deletion */
                                            Allele.create(String.valueOf(refBase), true)
                                            );
	                                
	                                v.add(new Base2(
                                            Allele.create(sb.toString().toUpperCase(),false),
                                            (byte)(qual/ce.getLength()),
                                            this.sample2index.get(sampleName),
                                            rec.getReadNegativeStrandFlag()
                                            ));
	                                }
                                readPos+=ce.getLength();
                                break;
                                }
                            case EQ: case M: case X:
                                {
                                for(int i=0; i< ce.getLength();++i)
                                    {
                                    MyVariantContext v=findContext(
                                            rec.getReferenceIndex(),
                                            refPos0 + i,
                                            Allele.create(String.valueOf(genomicSeq.charAt( refPos0 + i)), true)
                                            );
                                    

                                    v.add(new Base2(
                                            Allele.create(String.valueOf((char)bases[ readPos + i ]),false),
                                            quals[ readPos + i],
                                            this.sample2index.get(sampleName),
                                            rec.getReadNegativeStrandFlag()
                                            ));
                                    }
                                readPos+=ce.getLength();
                                refPos0+= ce.getLength();
                                break;
                                }

                            default : throw new
                            IllegalStateException("Case statement didn't deal with cigar op: "+ op);
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
            iter.close();
            this.variantContextWriter.close();
            return RETURN_OK;
            }
        catch (Exception e)
            {
        	LOG.error(e);
            return -1;
            }
        finally
            {
        	for(SamReader r:readers) CloserUtil.close(r);
            CloserUtil.close(this.indexedFastaSequenceFile);
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
