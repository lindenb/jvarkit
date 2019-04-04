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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.samtools.SamRecordJEXLFilter;
import htsjdk.variant.vcf.VCFIterator;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;



/**

BEGIN_DOC

## Motivation

Pour Sandro B.

GATK ClusteredReadPosition https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_ClusteredReadPosition.php  only works with Mutect2

The program looks for SNV in the VCF, go back to the reads in the bam.

For one variant , if all the reads contain the variant at less than 'distance' then the genotype is FILTERED

if all the reads are FILTERED, the variant is FILTERED


## Example

```
java -jar dist/vcfclusteredreadedge.jar -B in.bam in.vcf
```

```
find . -name "*.bam" > tmp.list 
java -jar dist/vcfclusteredreadedge.jar -B tmp.list  in.vcf
```


END_DOC
*/


@Program(name="vcfclusteredreadedge",
description="Variant annotation :   variants clustered near the ends of reads",
keywords={"sam","bam","vcf"}
)
public class VcfClusteredReadEdge extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfClusteredReadEdge.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-d","--distance"},description="minimal distance to the end of the **CLIPPED** read.")
	private int distance = 1 ;
	@Parameter(names={"-B","--bams"},description="path of indexed BAM path with read Groups. You can put those paths in a text file having a *.list sufffix")
	private List<String> bamList=new ArrayList<>();
	@Parameter(names={"-filter","--filter"},description=SamRecordJEXLFilter.FILTER_DESCRIPTION,converter=SamRecordJEXLFilter.StringConverter.class)
	private SamRecordFilter filter  = SamRecordJEXLFilter.buildDefault();
	@Parameter(names={"-gt","--gt"},description="Genotype FILTER name")
	private String genotypeFilterName="EDGEVAR";
	@Parameter(names={"-vt","--vt"},description="Variant FILTER name: set if ALL Genotypes have a variant near the edge.")
	private String variantFilterName="EDGEVAR";

	
	@Override
	protected int doVcfToVcf(final String inputName,final VCFIterator in,final VariantContextWriter out) {
		final List<File> bamFiles=  IOUtils.unrollFiles2018(this.bamList);
		final Map<String,List<SamReader>> sample2bam=new HashMap<>(bamFiles.size());
		final SamReaderFactory srf = super.createSamReaderFactory();
		try {
			final VCFHeader header=in.getHeader();
			for(final File bamFile: bamFiles)
				{
				LOG.info("Reading header for "+bamFile);
				final SamReader reader=srf.open(bamFile);
				if(!reader.hasIndex())
					{
					LOG.error("No BAM index available for "+bamFile);
					return -1;
					}		
				final SAMFileHeader samHeader=reader.getFileHeader();
				for(final SAMReadGroupRecord g:samHeader.getReadGroups())
					{
					if(g.getSample()==null) continue;
					final String sample=g.getSample();
					if(StringUtil.isBlank(sample)) continue;
					List<SamReader> readers = sample2bam.get(sample);
					if(readers==null)
						{
						readers=new ArrayList<>();
						sample2bam.put(sample,readers);
						}
					readers.add(reader);
					}
				}
			
			
			final VCFHeader h2=new VCFHeader(header);
			h2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));
			h2.addMetaDataLine(new VCFFilterHeaderLine(this.variantFilterName,
					"All genotypes have the variant at a distance of less or equal than "+this.distance+" bases in the reads."));

			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
			
			out.writeHeader(h2);
			while(in.hasNext())
				{
				final VariantContext ctx = progress.watch(in.next());
				
				final List<Genotype> genotypes = new ArrayList<>(ctx.getNSamples());
				int count_genotype_ok_edge=0;
				int count_genotype_bad_edge=0;
				for(int i=0;i< ctx.getNSamples();++i)
					{
					Genotype genotype = ctx.getGenotype(i);
					if(genotype.isHomRef() || genotype.isNoCall())
						{
						genotypes.add(genotype);
						continue;
						}
					
					final Set<Character> altBases= genotype.getAlleles().stream().
							filter(A->!(A.isReference() || A.isNoCall() || A.length()!=1 || A.isSymbolic())).
							map(A->A.getDisplayString().toUpperCase().charAt(0)).
							collect(Collectors.toSet())
							;
					if(altBases.isEmpty())
						{
						genotypes.add(genotype);
						continue;
						}
					boolean found_one_read_with_mutation=false;
					boolean all_positions_are_bad=true;
					final String sample = genotype.getSampleName();
					final List<SamReader> samReaders = sample2bam.get(sample);
					if(samReaders==null || samReaders.isEmpty())
						{
						genotypes.add(genotype);
						continue;
						}
					
					for(final SamReader sr: samReaders)
						{
						final SAMRecordIterator iter=sr.query(ctx.getContig(), ctx.getStart(), ctx.getEnd(), false);
						while(iter.hasNext())
							{
							final SAMRecord rec=iter.next();
							if(rec.getReadUnmappedFlag()) continue;
							if(filter.filterOut(rec)) continue;
							final SAMReadGroupRecord rg=rec.getReadGroup();
							if(!sample.equals(rg.getSample())) continue;
							final Cigar theCigar =rec.getCigar();
							if(theCigar==null) continue;
							int refPos= rec.getAlignmentStart();
							int readpos =0;
							String readBases=rec.getReadString();
							if(readBases.equals(SAMRecord.NULL_SEQUENCE_STRING)) continue;
							final List<CigarElement> cigars = new ArrayList<>(theCigar.getCigarElements());
							
							// remove 3' cigar clipping
							while(!cigars.isEmpty())
								{
								final CigarElement ce = cigars.get(cigars.size()-1);
								final CigarOperator op=ce.getOperator();
								if(op.equals(CigarOperator.SOFT_CLIP))
									{
									cigars.remove(cigars.size()-1);
									readBases=readBases.substring(0,readBases.length()-ce.getLength());
									}
								else if(op.equals(CigarOperator.HARD_CLIP))
									{
									cigars.remove(cigars.size()-1);
									}
								else
									{
									break;
									}
								}
							// remove 5' cigar clipping
							while(!cigars.isEmpty())
								{
								final CigarElement ce = cigars.get(0);
								final CigarOperator op=ce.getOperator();
								if(op.equals(CigarOperator.SOFT_CLIP))
									{
									cigars.remove(0);
									readBases=readBases.substring(ce.getLength());
									}
								else if(op.equals(CigarOperator.HARD_CLIP))
									{
									cigars.remove(0);
									}
								else
									{
									break;
									}
								}
							if(cigars.isEmpty()) continue;
							int position0_of_mitation_in_unclipped_read=-1;
							for(final CigarElement ce:cigars)
								{
								if( refPos > ctx.getEnd() ) break;
								final CigarOperator op=ce.getOperator();
								if(op.consumesReferenceBases() &&
									op.consumesReadBases() &&
									refPos>=ctx.getStart()
									)
									{
									for(int n=0;n< ce.getLength() && readpos+n < readBases.length() ;++n )
										{
										if(refPos+n!=ctx.getStart()) continue;
										final char readBase = readBases.charAt(readpos+n);
										if(altBases.contains(readBase))
											{
											found_one_read_with_mutation = true;
											int left_pos = readpos+n;
											int right_pos = readBases.length()-left_pos;
											position0_of_mitation_in_unclipped_read = Math.min(left_pos, right_pos);
											}
										break;
										}
									}
								if(op.consumesReferenceBases()) {
									refPos+= ce.getLength();
									}
								if(op.consumesReadBases()) {
									readpos+= ce.getLength();
									}
								}
							if(position0_of_mitation_in_unclipped_read==-1) continue;
							if(position0_of_mitation_in_unclipped_read +1 /* array 1 based */ <= this.distance) continue;
							all_positions_are_bad  = false;
							break;
							}
						iter.close();
						if(found_one_read_with_mutation && !all_positions_are_bad) break;
						}
					if(found_one_read_with_mutation && all_positions_are_bad)
						{
						genotype = new GenotypeBuilder(genotype).
							filter(this.genotypeFilterName).
							make();
						count_genotype_bad_edge++;
						}
					else
						{
						count_genotype_ok_edge++;
						}
						
					genotypes.add(genotype);
					}//end of for-each genotype
				
				// somethingWasChanged
				if(count_genotype_bad_edge>0)
					{
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					vcb.genotypes(genotypes);
					if(count_genotype_ok_edge==0) vcb.filter(this.variantFilterName);
					out.add(vcb.make());
					}
				else
					{
					out.add(ctx);
					}
				}
			progress.finish();
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{	
			sample2bam.values().stream().flatMap(L->L.stream()).forEach(R->CloserUtil.close(R));
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.distance< 0) {
			LOG.error("bad distance.");
			return -1;
		}
		
		if(StringUtil.isBlank(this.genotypeFilterName)) {
			LOG.error("bad genotypeFilterName.");
			return -1;
		}
		
		if(StringUtil.isBlank(this.variantFilterName)) {
			LOG.error("bad variantFilterName.");
			return -1;
		}
		try {
			return doVcfToVcf(args,this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	public static void main(final String[] args) {
		new VcfClusteredReadEdge().instanceMainWithExit(args);
	}

}
