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
package com.github.lindenb.jvarkit.tools.structvar;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.PeekIterator;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC


END_DOC

**/

@Program(name="saminversions",
	description="Scan inversions in SAM using paired read on same strand.",
	keywords={"sam","bam","sv","inversion"},
	modificationDate="20201211",
	creationDate="20200609",
	generate_doc=false
	)
public class SamInversions extends Launcher
	{
	private static final Logger LOG = Logger.build(SamInversions.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=CRAM_INDEXED_REFENCE)
	private Path referenceFaidx = null;
	@Parameter(names={"--max-size"},description="max distance between read and mate.",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int max_distance = 1_000 ;
	@Parameter(names={"-B","--bed","-r","--rgn"},description=IntervalListProvider.OPT_DESC,splitter= NoSplitter.class,converter=IntervalListProvider.StringConverter.class)
	private IntervalListProvider intervallistProvider = null;
	@Parameter(names={"-s","-supporting"},description="Don't print the variant if max FORMAT/DP < 's'")
	private int min_supporting_reads = 1;
	@Parameter(names={"--mapq"},description="min MAPQ")
	private int mapq = 30;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();


	
	@Override
	public int doWork(final List<String> args) {
		if(this.max_distance<100) {
			LOG.error("max size inversion must be >=100");
			return -1;
			}
	 	VariantContextWriter vcw= null;
		try {
			final Allele REF = Allele.create("N", true);
			final Allele INV = Allele.create("<INV>", false);
			

			final String input = super.oneFileOrNull(args);
			final Path inputPath = input==null?null:Paths.get(input);
			try(SamReader samReader= super.createSamReaderFactory().
					validationStringency(ValidationStringency.LENIENT).
					referenceSequence(this.referenceFaidx).open(input==null?
						SamInputResource.of(stdin()):
						SamInputResource.of(inputPath)
						)) {
				final SAMFileHeader samHeader = samReader.getFileHeader();
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(samHeader);

				final String sampleName= samHeader.getReadGroups().stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtil.isBlank(S)).
						findFirst().
						orElse(inputPath==null?"SAMPLE":IOUtils.getFilenameWithoutCommonSuffixes(inputPath));
				
				final ToIntFunction<SAMRecord> mateEnd = REC->SAMUtils.getMateCigar(REC)!=null?
						SAMUtils.getMateAlignmentEnd(REC):
						REC.getMateAlignmentStart();
				
				final ToIntFunction<SAMRecord> recordLen = REC->{
					final int start = Math.min(REC.getStart(), REC.getMateAlignmentStart());
					final int end = Math.max(REC.getEnd(), mateEnd.applyAsInt(REC));
					return CoordMath.getLength(start, end);
					};
				
				final QueryInterval queryIntervals[] = this.intervallistProvider==null?
						null:
						this.intervallistProvider.
							dictionary(dict).
							optimizedQueryIntervals()
							;
				
				final Set<VCFHeaderLine> meta=new HashSet<>();
				VCFStandardHeaderLines.addStandardFormatLines(meta,true,
						VCFConstants.GENOTYPE_KEY,
						VCFConstants.DEPTH_KEY
						
						);
				VCFStandardHeaderLines.addStandardInfoLines(meta,true,
						VCFConstants.DEPTH_KEY,
						VCFConstants.END_KEY
						);
				meta.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"SV length"));
				meta.add(new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String,"Structural variant type"));
				
				
				
				final VCFHeader header=new VCFHeader(meta,Collections.singleton(sampleName));
				JVarkitVersion.getInstance().addMetaData(this, header);
				header.setSequenceDictionary(dict);

				
				vcw = this.writingVariants.dictionary(dict).open(this.outputFile);
				vcw.writeHeader(header);

				final Predicate<SAMRecord> acceptRecord = rec->{
					if(!SAMRecordDefaultFilter.accept(rec, this.mapq)) return false;
					if(!rec.getReadPairedFlag()) return false;
					if(rec.getMateUnmappedFlag()) return false;
					if(!rec.getReferenceIndex().equals(rec.getMateReferenceIndex())) return false;
					if(rec.getReadNegativeStrandFlag()!=rec.getMateNegativeStrandFlag()) return false;
					final int len =  recordLen.applyAsInt(rec);
					if(len > max_distance)  return false;
					return true;
					};
				final List<Allele> alleles = Arrays.asList(REF,INV);
				final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.
						newInstance().
						dictionary(dict).
						logger(LOG).
						build();
				try(final CloseableIterator<SAMRecord> iter0 = (queryIntervals==null?samReader.iterator():samReader.query(queryIntervals,false))) {
					final PeekIterator<SAMRecord> iter = new PeekIterator<>(iter0);
					while(iter.hasNext()) {
						final SAMRecord rec = progress.apply(iter.next());
						if(!acceptRecord.test(rec)) continue;
						final int invStart = rec.getAlignmentStart();
						int invEnd = Math.max(rec.getEnd(),mateEnd.applyAsInt(rec));
						if(CoordMath.getLength(invStart, invEnd)>this.max_distance) continue;
						int dp = 1;
						while(iter.hasNext()) {
							final SAMRecord rec2 = iter.peek();
							if(!acceptRecord.test(rec2)) {
								iter.next();
								continue;
								}
							if(rec2.getAlignmentStart() > invEnd ||
								!rec2.contigsMatch(rec) ||
								CoordMath.getLength(invStart, rec2.getAlignmentStart())>this.max_distance) {
								break;
								}
							iter.next();//consumme
							invEnd = Math.max(invEnd,Math.max(rec2.getEnd(),mateEnd.applyAsInt(rec2)));
							dp++;
							}
						
						if(dp >= this.min_supporting_reads) {
							final VariantContextBuilder vcb = new VariantContextBuilder();
							
							
							vcb.chr(rec.getContig());
							vcb.start(invStart);
							vcb.stop(invEnd);
							vcb.alleles(alleles);
							vcb.attribute(VCFConstants.END_KEY, invEnd);
							
							final GenotypeBuilder gb = new GenotypeBuilder(sampleName,alleles);
							gb.DP(dp);					
								
							
							vcb.genotypes(Collections.singletonList(gb.make()));
							vcb.alleles(alleles);
							vcb.attribute(VCFConstants.DEPTH_KEY, dp);
							vcb.attribute(VCFConstants.SVTYPE, "INV");
							vcb.attribute("SVLEN", CoordMath.getLength(invStart,invEnd));
							vcw.add(vcb.make());
							}
						}
					}
				progress.close();
				vcw.close();
				}
			return 0;
			} 
		catch (final Throwable e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcw);	
			}
		}

	
	public static void main(final String[] args)
		{
		new SamInversions().instanceMainWithExit(args);
		}
	}
