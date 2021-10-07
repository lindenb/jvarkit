/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.structvar.indexcov;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Motivation

Same as indexcov2sv but merge contigous segments.

## Input

input is a tab-delimited file created by e.g: indexcov (https://github.com/brentp/goleft/tree/master/indexcov)

```
#chrom  start  end     SampleBB  SampleBC  SampleBD  SampleBE  SampleBF  SampleBG  SampleBH
chr1    23778  40778   1.59      1.31      1.67      1.61      1.83      1.52      1.48
chr1    29106  46106   1.9       1.54      1.72      1.97      1.88      1.53      1.95
chr1    84581  101581  0.764     0.841     1.2       1.16      1.18      1.13      1.23
chr1    15220  32220   0.355     0.704     1.09      0.784     0.81      1.37      0.954
chr1    58553  75553   0.353     0.436     0.912     0.836     1.16      1.09      0.611
chr1    19347  36347   0.381     0.411     0.811     0.795     1.16      1.22      0.495
chr1    81062  98062   1.09      0.972     1.35      1.22      1.66      1.76      1.1
chr1    17353  34353   1.06      1.06      1.23      1.26      1.44      1.43      1.03
chr1    48498  65498   1.08      0.996     1.28      1.44      1.52      1.57      1.05
```


END_DOC
 */
@Program(
		name="indexcov2sv",
		description="same as indexcov2vcf but merge segments",
		keywords={"cnv","duplication","deletion","sv","indexcov"},
		creationDate="20200228",
		modificationDate="20211007"
		)
public class IndexCovToSV extends Launcher {
	private static final Logger LOG = Logger.build(IndexCovToSV.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-t","--treshold"},description=IndexCovUtils.TRESHOLD_OPT_DESC+". " +FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double treshold = IndexCovUtils.DEFAULT_TRESHOLD;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refFile = null;
	@Parameter(names={"-d","--distance"},description="Merge segment if they are within that distance. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int merge_distance=0;
	@Parameter(names={"-norm","--norm"},description="normalize: only one ALT allele per variant")
	private boolean normalize_one_allele = false;

	
	
	@ParametersDelegate
	private WritingSortingCollection sortingDelegate = new WritingSortingCollection();
	@ParametersDelegate
	private WritingVariantsDelegate writingDelegate = new WritingVariantsDelegate();

	
	private static class AltAllele {
		final IndexCovUtils.SvType type;
		final int end;
		final Allele allele;
		int ac = 0;
		AltAllele(IndexCovUtils.SvType type,int end,final Allele allele) {
			this.type = type;
			this.end= end;
			this.allele = allele;
			}
		
		}
	
	private static class Cell {
		int tid;
		int start;
		int end;
		int sample_id;
		double score;
		
		int sortChromStart(final Cell o) {
			int i = Integer.compare(this.tid,o.tid);
			if(i!=0) return i;
			i = Integer.compare(this.start,o.start);
			return i;
			}
		
		int sortSamplePosition(final Cell o) {
			final int i = Integer.compare(this.sample_id,o.sample_id);
			if(i!=0) return i;
			return sortChromStart(o);
			}
		
		int sortPositionSample(final Cell o) {
			final int i = sortChromStart(o);
			if(i!=0) return i;
			return Integer.compare(this.sample_id,o.sample_id);
			}
		}
	
	private static class CellCodec extends AbstractDataCodec<Cell> {
		@Override
		public Cell decode(DataInputStream dis) throws IOException {
			final Cell c= new Cell();
			try {
				c.tid = dis.readInt();
				c.start = dis.readInt();
				c.end = dis.readInt();
				c.sample_id = dis.readInt();
				c.score  = dis.readDouble();
				} catch(final EOFException err) {
					return null;
				}
			return c;
			}
		@Override
		public void encode(DataOutputStream out, Cell o) throws IOException {
			out.writeInt(o.tid);
			out.writeInt(o.start);
			out.writeInt(o.end);
			out.writeInt(o.sample_id);
			out.writeDouble(o.score);
			}
		
		@Override
		public AbstractDataCodec<Cell> clone() {
			return new CellCodec();
		}
	}
	
	public IndexCovToSV() {
	}
	

	
	
	@Override
	public int doWork(final List<String> args) {
		final IndexCovUtils indexCovUtils  = new IndexCovUtils(this.treshold);
		final CharSplitter tab = CharSplitter.TAB;
		
		
		try {
			
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refFile);
			
			final List<String> samples;
			SortingCollection<Cell> sorter1 = SortingCollection.newInstance(
					Cell.class,
					new CellCodec(),
					(A, B)->A.sortSamplePosition(B),
					this.sortingDelegate.getMaxRecordsInRam(),
					this.sortingDelegate.getTmpPaths()
					);
			sorter1.setDestructiveIteration(true);

			
			try(BufferedReader r =  super.openBufferedReader(oneFileOrNull(args))) {
				String line = r.readLine();
				if(line==null) {		
					LOG.error( "Cannot read first line of input");
					return -1;
					}
				String tokens[] = tab.split(line);
				if(tokens.length<4 ||
					!tokens[0].equals("#chrom") ||
					!tokens[1].equals("start") ||
					!tokens[2].equals("end")) {
					LOG.error( "bad first line "+line );
					return -1;
					}
				samples = Arrays.asList(tokens). subList(3,tokens.length);			
					
					
				while((line=r.readLine())!=null) {
					if(StringUtil.isBlank(line)) continue;
					tokens =  tab.split(line);
					if(tokens.length!=3+samples.size()) {
						throw new JvarkitException.TokenErrors(samples.size()+3, tokens);
						}
					final SAMSequenceRecord ssr = dict.getSequence(tokens[0]);
					if(ssr==null) {
						LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(tokens[0],dict));
						return -1;
						}
					for(int i=3;i< tokens.length;i++) {
						final double score = Double.parseDouble(tokens[i]);
						if(indexCovUtils.getType(score).isAmbigous()) continue;
						final Cell cell = new Cell();
						cell.tid = ssr.getSequenceIndex();
						cell.start = Integer.parseInt(tokens[1]);
						cell.end = Integer.parseInt(tokens[2]);
						cell.sample_id = i-3;
						cell.score = score;
						sorter1.add(cell);
						}
					}// while readline
					
					sorter1.doneAdding();
				}
				/* merge adjacent blocks */
				
				SortingCollection<Cell> sorter2 = SortingCollection.newInstance(
						Cell.class,
						new CellCodec(),
						(A, B)->A.sortPositionSample(B),
						this.sortingDelegate.getMaxRecordsInRam(),
						this.sortingDelegate.getTmpPaths()
						);
				sorter2.setDestructiveIteration(true);
				
				
				try(CloseableIterator<Cell> iter1 = sorter1.iterator()) {
					try(PeekableIterator<Cell> peeker1 = new PeekableIterator<>(iter1)) {
						while(peeker1.hasNext()) {
							final Cell first = peeker1.next();
							while(peeker1.hasNext()) {
								final Cell second = peeker1.peek();
								if(first.sample_id!=second.sample_id) break;
								if(first.tid!=second.tid) break;
								if(first.end+(this.merge_distance+1)<second.start) break;
								if(!indexCovUtils.getType(first.score).equals(indexCovUtils.getType(second.score))) break;
								first.end = second.end; // extend first with end of second
								peeker1.next();//consumme
								}
							sorter2.add(first);
							}// while peeker1
						}
					}
				sorter1.cleanup();sorter1=null;
				sorter2.doneAdding();
				
				final Set<VCFHeaderLine> metaData = new HashSet<>();
				VCFStandardHeaderLines.addStandardFormatLines(metaData, true,
						VCFConstants.GENOTYPE_KEY,
						VCFConstants.GENOTYPE_QUALITY_KEY
						);
				VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.END_KEY,
						VCFConstants.ALLELE_FREQUENCY_KEY,
						VCFConstants.ALLELE_COUNT_KEY,
						VCFConstants.ALLELE_NUMBER_KEY
						);
				
				final VCFFilterHeaderLine filterHomDel = new VCFFilterHeaderLine("HOM_DEL", "There is one Homozygous deletion.");
				metaData.add(filterHomDel);
				final VCFFilterHeaderLine filterHomDup = new VCFFilterHeaderLine("HOM_DUP", "There is one Homozygous duplication.");
				metaData.add(filterHomDup);

				/** raw value in indexcov */
				final VCFFormatHeaderLine foldHeader = new VCFFormatHeaderLine("F", 1, VCFHeaderLineType.Float,"Relative number of copy: 0.5 deletion 1 normal 2.0 duplication");
				metaData.add(foldHeader);
				final VCFInfoHeaderLine infoSvType = new VCFInfoHeaderLine(VCFConstants.SVTYPE, 1, VCFHeaderLineType.String,"SV type");
				metaData.add(infoSvType);

				final VCFInfoHeaderLine infoSvLen = new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"ABS(SV length)");
				metaData.add(infoSvLen);

				/*
				final VCFInfoHeaderLine infoSvLens = new VCFInfoHeaderLine("SVLENS", 1, VCFHeaderLineType.Integer,"Allele length");
				metaData.add(infoSvLens);
				final VCFInfoHeaderLine infoSvEnds = new VCFInfoHeaderLine("SVENDS", 1, VCFHeaderLineType.Integer,"Allele ends");
				metaData.add(infoSvEnds);
				 */

				
				
				final VCFHeader vcfHeader = new VCFHeader(metaData, samples);
		  		JVarkitVersion.getInstance().addMetaData(this, vcfHeader);
				vcfHeader.setSequenceDictionary(dict);
				try(VariantContextWriter vcw = this.writingDelegate.dictionary(dict).open(outputFile)) {
					vcw.writeHeader(vcfHeader);
	
					final Allele REF_ALLELE =Allele.create("N",true);
					
					try(CloseableIterator<Cell> iter2 = sorter2.iterator()) {
						try(EqualRangeIterator<Cell> peeker2 = new EqualRangeIterator<>(iter2,(A,B)->A.sortChromStart(B))) {
						while(peeker2.hasNext()) {
							final List<Cell> array0 = peeker2.next();
							final Set<Integer> distinct_ends = array0.stream().map(C->C.end).collect(Collectors.toCollection(TreeSet::new));
							final List<List<Cell>> array_of_array;
							if(this.normalize_one_allele && distinct_ends.size()>1) {
								array_of_array = new ArrayList<>(distinct_ends.size());
								for(final int end_variant: distinct_ends) {
									array_of_array.add(array0.stream().filter(C->C.end==end_variant).collect(Collectors.toList()));
									}
								}
							else
								{
								array_of_array = Collections.singletonList(array0);
								}
							
							for(final List<Cell> array: array_of_array) {
								int an=0;
								if(array.stream().
										map(C->indexCovUtils.getType(C.score)).
										noneMatch(C->C.isVariant())
										) continue;
								
								final Cell first = array.get(0);
								final int max_end = array.stream().mapToInt(C->C.end).max().getAsInt();
								final List<Genotype> genotypes = new ArrayList<>(array.size());
								final VariantContextBuilder vcb = new VariantContextBuilder();
								vcb.chr(dict.getSequence(first.tid).getSequenceName());
								vcb.start(first.start);
								vcb.stop(max_end);
								vcb.attribute(VCFConstants.END_KEY, max_end);
			
								final List<AltAllele> altAlleles = new ArrayList<>();
								boolean has_hom_del = false;
								boolean has_hom_dup = false;

								
								for(int i=0;i< array.size();i++) {
									final Cell cell = array.get(i);
									final String sampleName = samples.get(cell.sample_id);
									an+=2;
									final GenotypeBuilder gb;
									final IndexCovUtils.SvType type = indexCovUtils.getType(cell.score);
									final Allele allele;
									AltAllele altAllele;
									if(type.isDeletion() || type.isDuplication()) {
										altAllele = altAlleles.stream().filter(C->C.end==cell.end && C.type.equals(type)).findFirst().orElse(null);
										if(altAllele==null) {
											allele = Allele.create("<"+(type.isDeletion()?"DEL":"DUP")+"."+(1+altAlleles.size())+">", false);
											altAllele = new AltAllele(type,cell.end,allele);
											altAlleles.add(altAllele);
											}
										else
											{
											allele = altAllele.allele;
											}
										}
									else
										{
										altAllele = null;
										allele = null;
										}
									
									
									
									switch(type) {			
										case REF :
											gb = new GenotypeBuilder(sampleName,Arrays.asList(REF_ALLELE,REF_ALLELE));
											break;
										case HET_DEL:
											gb = new GenotypeBuilder(sampleName,Arrays.asList(REF_ALLELE,allele));
											altAllele.ac++;
											break;
										case HOM_DEL:
											gb = new GenotypeBuilder(sampleName,Arrays.asList(allele,allele));						
											vcb.filter(filterHomDel.getID());
											altAllele.ac+=2;
											has_hom_del = true;
											break;
										case HET_DUP:
											gb = new GenotypeBuilder(sampleName,Arrays.asList(REF_ALLELE,allele));
											altAllele.ac++;
											break;
										case HOM_DUP:
											gb = new GenotypeBuilder(sampleName,Arrays.asList(allele,allele));
											vcb.filter(filterHomDup.getID());
											altAllele.ac+=2;
											has_hom_dup = true;
											break;
										default:
											gb = new GenotypeBuilder(sampleName,Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
											break;
										}
									
									gb.attribute(foldHeader.getID(),cell.score);
									
									if(!type.isAmbigous()) {
										gb.GQ(type.getGenotypeQuality(cell.score));
										}
														
									genotypes.add(gb.make());	
									}//end loop over genotypes
								
								
								
								final List<Allele> alleles = altAlleles.stream().map(A->A.allele).collect(Collectors.toCollection(ArrayList::new));
								if(has_hom_del) {
									vcb.filter(filterHomDel.getID());
								} else if(has_hom_dup) {
									vcb.filter(filterHomDup.getID());
								} else
								{
									vcb.passFilters();
								}
								
								alleles.add(0,REF_ALLELE);
								vcb.alleles(alleles);
								vcb.genotypes(genotypes);
								if(altAlleles.stream().allMatch(C->C.type.isDeletion())) {
									vcb.attribute(VCFConstants.SVTYPE,"DEL");
								} else if(altAlleles.stream().allMatch(C->C.type.isDuplication())) {
									vcb.attribute(VCFConstants.SVTYPE,"DUP");
								} else  {
									vcb.attribute(VCFConstants.SVTYPE,"MIXED");
								} 
								
								vcb.attribute(infoSvLen.getID(), altAlleles.stream().mapToInt(C->1+C.end-first.start).max().orElse(0));
								vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,altAlleles.stream().map(A->A.ac).collect(Collectors.toList()));
								vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,an);
								if(an>0) {
									final int final_an= an;
									vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,
											altAlleles.stream().map(A->A.ac/(double)final_an).collect(Collectors.toList())
											);
									}
								
								vcw.add(vcb.make());
								}
							} //end loop while iter
						}
					}
				sorter2.cleanup();
				sorter2=null;
				}
			
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}

		public static void main(final String[] args) {
			new IndexCovToSV().instanceMainWithExit(args);
			}

}
