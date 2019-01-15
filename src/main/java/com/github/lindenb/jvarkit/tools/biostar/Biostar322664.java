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


*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFFileReader;

/** 

BEGIN_DOC

## Input

  * BAM : MUST be sorted using Picard SortSam (see https://github.com/samtools/hts-specs/issues/5 )
  * VCF : only SNP are considered. Genotypes are ignored (all ALT alleles are observed regardless of the sample/genotype )

##Example

```
$ java -jar picard.jar SortSam I=src/test/resources/S1.bam O=query.bam SO=queryname
$ java -jar dist/biostar322664.jar -V src/test/resources/S1.vcf.gz query.bam  


RF02_358_926_2:0:0_2:1:0_83	83	RF02	857	60	70M	=	358	-569	GACGTGAACTATATAATTAAAATGGACAGAAATCTGCCATCAACAGCTAGATATATAAGACCTAATTTAC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
RF02_362_917_2:0:0_2:1:0_6f	147	RF02	848	60	70M	=	362	-556	ATAAGGAATCACGTTAACTATATACTTAAAATGGACTGAAATCTGCCATCAACAGCTAGATATATAAGAC	2222222222222222222222222222222222222222222222222222222222222222222222	RG:Z:S1	NM:i:3	AS:i:55	XS:i:0
(...)
```

```
$ java -jar dist/biostar322664.jar -nm  -index -V src/test/resources/S1.vcf.gz src/test/resources/S1.bam 
(...)
RF02_827_1385_4:1:0_2:0:0_3f    163     RF02    827     60      70M     =       1316    559     ATCAATTACATTCCTGAAAGGATAAGGAATGAGGTTAACTATCTACTTAAAATGGACAGAAATCTGCCAA  2222222222222222222222222222222222222222222222222222222222222222222222  RG:Z:S1 NM:i:5  AS:i:53XS:i:0
RF02_827_1292_4:0:0_1:0:0_50    163     RF02    827     60      70M     =       1223    466     TTCAATTACATTCCTGCAAGGATAAGGAATGCCGTTAACTATATACTTAATAAGGACAGAAATCTGGCAT  2222222222222222222222222222222222222222222222222222222222222222222222  RG:Z:S1 NM:i:4  AS:i:51XS:i:0
(...)
```


END_DOC

*/

@Program(name="biostar322664",
	description="Extract PE Reads (with their mates) supporting variants in vcf file",
	keywords= {"sam","bam","vcf"},
	biostars=322664
	)
public class Biostar322664 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar322664.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-V","--variant"},description="Variant VCF file. This tool **doesn't work** with INDEL/SV.",required=true)
	private File vcfFile = null;
	@Parameter(names={"-nm","--no-mate"},description="Disable the 'mate' function. BAM is not expected to be sorted with picard (bam can be sorted on coordinate), but mate will not be written.")
	private boolean input_is_not_sorted_on_queryname = false;
	@Parameter(names={"-index","--index"},description="Use the VCF input to query the BAM using bai index. Faster for large bam + small VCF. Require option `--no-mate` and the bam file to be indexed. ")
	private boolean use_bam_index = false;
	@Parameter(names={"-x","-X"},description="If defined, add the variant(s) information in a 'X' metadata. One character only.")
	private String meta_tag_str = null;
	@Parameter(names={"-all","--all"},description="used in complement of option -X . Will output any SamRecord, but some reads will carrying the information about the variant in a 'Xx' attribute.")
	private boolean output_all = false;
	@Parameter(names={"-pair","--pair"},description="pair mode: the paired read and it's pair muts BOTH carry at least one variant")
	private boolean both_pair_mode = false;

	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();

	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.use_bam_index && !this.input_is_not_sorted_on_queryname)
			{
			LOG.error("Cannot use option -index without option -nm");
			return -1;
			}
		if(this.both_pair_mode && this.input_is_not_sorted_on_queryname) {
			LOG.error("Cannot use option -and without sorting on queryname");
			return -1;
			}
		/* will be used as a SAMRecord attribute 'X'+meta_char containing the data about the variants */
		final char meta_char;
		if(this.meta_tag_str!=null) {
			if(this.meta_tag_str.length()!=1) {
				LOG.error("Meta tag should be only one character, got "+meta_tag_str);
				return -1;			
				}
			meta_char = this.meta_tag_str.charAt(0);
		} else
			{
			meta_char='\0';
			}
		
		if(this.output_all && meta_char=='\0') {
			LOG.error("Cannot output all if not meta char -X defined");
			return -1;
		}
		
		final IntervalTreeMap<VariantContext> variantMap=new IntervalTreeMap<>();
		VCFFileReader vcfFileReader=null;
		CloseableIterator<VariantContext> viter=null;
		SamReader samReader = null;
		SAMFileWriter samFileWriter = null;
		try {
			LOG.info("reading VCF "+this.vcfFile+" in memory...");
			vcfFileReader = new VCFFileReader(this.vcfFile, false);
			viter= vcfFileReader.iterator();
			viter.stream().
				filter(V->V.isVariant() && V.getReference().length()==1 && V.getAlternateAlleles().stream().anyMatch(A->A.length()==1)).
				map(V->new VariantContextBuilder(V).attributes(Collections.emptyMap()).noID().unfiltered().noGenotypes().make()).
				forEach(V->variantMap.put(new Interval(V.getContig(), V.getStart(), V.getEnd()),V));
			CloserUtil.close(viter);viter=null;
			vcfFileReader.close();vcfFileReader=null;
			LOG.info("Done Reading: "+this.vcfFile);
			

			samReader =super.openSamReader(oneFileOrNull(args));
			final SAMFileHeader header = samReader.getFileHeader();
			if(!this.input_is_not_sorted_on_queryname &&
				header.getSortOrder()!=SAMFileHeader.SortOrder.queryname) {
				LOG.error("Expected SAM input to be sorted on "+SAMFileHeader.SortOrder.queryname+" but got "+header.getSortOrder() +
						" See the options to work without 'mate'.");
				return -1;
				}
			samFileWriter = this.writingBamArgs.openSAMFileWriter(this.outputFile, header, true);
			final CloseableIterator<SAMRecord> iter ;
			if(!use_bam_index) {
				iter = samReader.iterator();
				}
			else
				{
				final SAMSequenceDictionary dict = header.getSequenceDictionary();
				if(dict==null) throw new JvarkitException.BamDictionaryMissing("input");
				final List<QueryInterval> intervals = new ArrayList<>();
				for(final VariantContext ctx:variantMap.values())
					{
					int tid = dict.getSequenceIndex(ctx.getContig());
					if(tid<0)  throw new JvarkitException.ContigNotFoundInDictionary(ctx.getContig(), dict);
					intervals.add(new QueryInterval(tid, ctx.getStart(), ctx.getEnd()));
					}
				QueryInterval intervalArray[] = intervals.toArray(new QueryInterval[intervals.size()]);
				intervalArray = QueryInterval.optimizeIntervals(intervalArray);
				iter = samReader.query(intervalArray, false);
				}
			
			
			Iterator<List<SAMRecord>> eq_range = null;
			if(this.input_is_not_sorted_on_queryname)
				{
				eq_range = new AbstractIterator<List<SAMRecord>>()
					{
					@Override
					protected List<SAMRecord> advance()
						{
						return iter.hasNext()?Collections.singletonList(iter.next()):null;
						}
					};
				}
			else
				{
				eq_range = new EqualRangeIterator<>(iter,(R1,R2)->R1.getReadName().compareTo(R2.getReadName()));
				}
			while(eq_range.hasNext()) {
				final List<SAMRecord> array = eq_range.next();
				final boolean carry_variant[]=new boolean[array.size()];
				Arrays.fill(carry_variant, false);
				
				for(int record_index=0;record_index < array.size();++record_index) {
					final SAMRecord rec= array.get(record_index);
					//final boolean debug = rec.getReadName().equals("SRR5229653.1238236");
					//if(debug==true) LOG.info("got read! "+rec+" array.size="+array.size());
					if(rec.getReadUnmappedFlag()) {
						continue;
					}
					
					final Set<String> meta_attribute = (meta_char=='\0'?null:new HashSet<>());
					final Collection<VariantContext> variants = variantMap.getOverlapping(rec);
					if(variants.isEmpty()) {
						continue;
					}
					
					final Cigar cigar=rec.getCigar();
					if(cigar==null || cigar.isEmpty()) continue;
					
					final byte bases[] = rec.getReadBases();
					if(bases==null || bases.length==0 || bases==SAMRecord.NULL_SEQUENCE) continue;
					
					for(final VariantContext ctx:variants) {
						if(ctx.getReference().length()!=1) continue;
						for(final Allele alt:ctx.getAlternateAlleles()) {
							if(alt.isNoCall() || alt.isSymbolic()) continue;
							if(alt.length()!=1) continue;
							final char altbase= (char)Character.toUpperCase(alt.getBases()[0]);
							int ref1 =rec.getUnclippedStart();
							int readpos = 0;
							
							for(final CigarElement ce:cigar) {
								if(ref1> ctx.getEnd()) {
									break;
								}
								final CigarOperator op =ce.getOperator();
								switch(op) {
								case P:break;
								case S: 
									{
									readpos+=ce.getLength();
									ref1+=ce.getLength();
									break;
									}
								case H:
									{
									ref1+=ce.getLength();
									break;
									}
								case I: readpos+=ce.getLength();break;
								case D:case N: ref1+=ce.getLength();break;
								case M:case EQ: case X: 
									{
									for(int x=0;x< ce.getLength();++x)
										{
										final int rp2 = readpos+x;
										if(rp2>=0 && rp2<bases.length && ref1+x==ctx.getStart())
											{
											final char readbase = Character.toUpperCase((char)bases[rp2]);
											
											if(readbase==altbase)
												{
												//if(debug) LOG.debug(rec.getReadName()+" read["+rp2+"]="+readbase+" ref["+ctx.getStart()+"]="+altbase);
												carry_variant[record_index]=true;
												if(meta_attribute!=null) {
													meta_attribute.add(ctx.getContig()+"|"+ctx.getStart()+"|"+ctx.getReference().getDisplayString()+"|"+ctx.getAlternateAlleles().stream().filter(A->A.length()==1).map(A->A.getDisplayString()).collect(Collectors.joining(",")));
													}
												break;
												}
											}
										}
									readpos+=ce.getLength();
									ref1+=ce.getLength();
									break;
									}
								default: throw new IllegalStateException("bad cigar operator "+op);
								}
							}
						}//end of cigar loop					
					if(meta_attribute!=null && !meta_attribute.isEmpty()) {
						rec.setAttribute("X"+meta_char, String.join(";",meta_attribute));
						}
					} // end of loop over variants
				}// end of for(Samrecord in array)
				
				boolean any_match=false;
				for(boolean m: carry_variant) {
					if(m) {
						any_match=true;
						break;
					}
				}
				
				if(any_match && this.both_pair_mode)  {
					any_match = false;//reset
					for(int x=0;x+1< array.size() && !any_match;++x) {
						if(!carry_variant[x]) continue;
						final SAMRecord r1 = array.get(x);
						if(!r1.getReadPairedFlag()) continue;
						for(int y=x+1;y< array.size() && !any_match;++y) {
							if(!carry_variant[y]) continue;
							// check r2 is mate of r1
							final SAMRecord r2 = array.get(x);
							if(!r2.getReadPairedFlag()) continue;
							if(r1.getFirstOfPairFlag()&& r2.getFirstOfPairFlag()) continue;
							if(!r1.getFirstOfPairFlag()&& !r2.getFirstOfPairFlag()) continue;
							if(r1.getAlignmentStart()!=r2.getMateAlignmentStart()) continue;
							if(r2.getAlignmentStart()!=r1.getMateAlignmentStart()) continue;
							if(!r1.getReferenceName().equals(r2.getMateReferenceName())) continue;
							any_match = true;
						}
					}
				}
				
				
				
				if(any_match || this.output_all) {
					for(final SAMRecord rec:array) samFileWriter.addAlignment(rec);
					}
				}
			CloserUtil.close(eq_range);
 			iter.close();
			samFileWriter.close();samFileWriter=null;
			samReader.close();samReader=null;
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
		} finally
			{
			CloserUtil.close(samReader);
			CloserUtil.close(samFileWriter);
			CloserUtil.close(vcfFileReader);
			}
	}
		
	
	public static void main(final String[] args) throws IOException
		{
		new Biostar322664().instanceMainWithExit(args);
		}
		

	}
