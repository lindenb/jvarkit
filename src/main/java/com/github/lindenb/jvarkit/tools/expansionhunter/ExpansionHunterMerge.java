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
package com.github.lindenb.jvarkit.tools.expansionhunter;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.variant.vcf.VCFRecordCodec;

/**
BEGIN_DOC
 
# Input
 
Input is a list of indexed vcf files or one file with the '.list' suffix containing the path to the vcfs
 

END_DOC
*/
@Program(name="expansionhuntermerge",
description="Merge Vcf from ExpansionHunter.",
keywords= {"vcf","merge","ExpansionHunter"},
creationDate="20210210",
modificationDate="20210210"
)
public class ExpansionHunterMerge extends Launcher {
	private static final Logger LOG = Logger.build( ExpansionHunterMerge.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();

	@Override
		public int doWork(List<String> args) {
		 	SortingCollection<VariantContext> sorter = null;
			try {
				final List<Path> inputs = IOUtils.unrollPaths(args);
				if(inputs.isEmpty()) {
					LOG.info("no input file.");
					return -1;
					}
				final String REPID="REPID";
				final String RU="RU";
				final String fakeSample="___FAKE";
				SAMSequenceDictionary dict = null;
				final Set<String> samples = new TreeSet<>();
				final Set<VCFHeaderLine> metaData = new HashSet<>();
				for(final Path path:inputs) {
					try(VCFReader r=VCFReaderFactory.makeDefault().open(path, false)) {
						final VCFHeader header = r.getHeader();
						if(header.getInfoHeaderLine(REPID)==null) {
							LOG.error("missing INFO/"+REPID);
							return -1;
							}
						final SAMSequenceDictionary d = header.getSequenceDictionary();
						if(d!=null) {
							if(dict==null) {
								dict = d;
								}
							else
								{
								SequenceUtil.assertSequenceDictionariesEqual(d, dict);
								}
							}
						if(header.getNGenotypeSamples()!=1) {
							LOG.error("expected one and only one genotyped sample in "+path);
							return -1;
							}
						final String sn = header.getSampleNamesInOrder().get(0);
						if(sn.equals(fakeSample)) {
							LOG.error(sn+" cannot be named "+fakeSample+" in "+path);
							return -1;
							}
						if(samples.contains(sn)) {
							LOG.error("duplicate sample "+sn+" in "+path);
							return -1;
							}
						metaData.addAll(header.getMetaDataInInputOrder());
						samples.add(sn);
						}
					}
				final VCFHeader tmpHeader = new VCFHeader(metaData, Arrays.asList(fakeSample));
				final VCFInfoHeaderLine sampleInfo = new VCFInfoHeaderLine("SRCSAMPLE",1,VCFHeaderLineType.String,"SRC SAMPLE");
				final Comparator<String> ctgComparator= dict==null?
						(A,B)->A.compareTo(B):
						new ContigDictComparator(dict)
						;
				
				final BiFunction<VariantContext, String, String> getAtt = (V,A)->{
					final String s1 = V.getAttributeAsString(A, "");
					if(StringUtils.isBlank(s1)) throw new IllegalStateException("INFO/"+A+" missing in "+V);
					return s1;
					};
				
				final Comparator<VariantContext> comparator1= (V1,V2) ->{
					String s1 = V1.getContig();
					String s2 = V2.getContig();
					int i= ctgComparator.compare(s1,s2);
					if(i!=0) return i;
				
					i= Integer.compare(V1.getStart(), V2.getStart());
					if(i!=0) return i;
					
					s1 = getAtt.apply(V1,REPID);
					s2 = getAtt.apply(V2,REPID);
					return s1.compareTo(s2);
					};
					
				final Comparator<VariantContext> comparator2 = (V1,V2) ->{
					final int i= comparator1.compare(V1, V2);
					if(i!=0) return i;
					final String s1 = getAtt.apply(V1,sampleInfo.getID());
					final String s2 = getAtt.apply(V2,sampleInfo.getID());
					return s1.compareTo(s2);
					};
				
				
				tmpHeader.addMetaDataLine(sampleInfo);
				sorter = SortingCollection.newInstance(
                        VariantContext.class,
                        new VCFRecordCodec(tmpHeader, false),
                        comparator2,
                        this.writingSortingCollection.getMaxRecordsInRam(),
                        this.writingSortingCollection.getTmpPaths()
                        );
				sorter.setDestructiveIteration(true);
				for(final Path path:inputs) {
					LOG.info("adding "+path);
					try(VCFReader r=VCFReaderFactory.makeDefault().open(path, false)) {
						try(CloseableIterator<VariantContext> iter = r.iterator()) {
							while(iter.hasNext()) {
								final VariantContext ctx = iter.next();
								
								final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
								final Genotype gt0 = ctx.getGenotype(0);
								final Genotype gt = new GenotypeBuilder(gt0).name(fakeSample).make();
								vcb.genotypes(Collections.singletonList(gt));
								vcb.attribute(sampleInfo.getID(), gt0.getSampleName());
								sorter.add(vcb.make());
								}
							}
						}
					}
				
				
				sorter.doneAdding();
				final VCFHeader outputHeader = new VCFHeader(metaData, samples);
				if(dict!=null) outputHeader.setSequenceDictionary(dict);
				JVarkitVersion.getInstance().addMetaData(this, outputHeader);
				try(VariantContextWriter w= writingVariants.dictionary(dict).open(this.outputFile)) {
					w.writeHeader(outputHeader);
					try(CloseableIterator<VariantContext> iter = sorter.iterator()) {
						final EqualRangeIterator<VariantContext> eq = new EqualRangeIterator<>(iter,comparator1);
						while(eq.hasNext()) {
							final List<VariantContext> calls = eq.next();
							final VariantContext first = calls.get(0);
							final Set<Allele> altAllelesSet = calls.stream().
									flatMap(V->V.getGenotypes().stream()).
									flatMap(GT->GT.getAlleles().stream()).
									filter(A->!(A.isReference() || A.isNoCall())).
									collect(Collectors.toSet());
								if(altAllelesSet.isEmpty()) continue;
								final List<Allele> altAllelesList = new ArrayList<>(altAllelesSet);
								final List<Allele> vcAlleles = new ArrayList<>(altAllelesList.size()+1);
								vcAlleles.add(first.getReference());
								vcAlleles.addAll(altAllelesList);

								
								final List<Genotype> genotypes = new ArrayList<>(samples.size());
								for(final String sn : samples) {
									final VariantContext vcs = calls.stream().
											filter(V->getAtt.apply(V,sampleInfo.getID()).equals(sn)).
											findFirst().orElse(null);
									final Genotype gt;
									if(vcs==null) {
										gt=GenotypeBuilder.createMissing(sn,2);
										}
									else
										{
										gt= new GenotypeBuilder(vcs.getGenotype(0)).
												name(sn).
												make();
										}
									genotypes.add(gt);
									}
								final VariantContextBuilder vcb = new VariantContextBuilder(null,
										first.getContig(),
										first.getStart(), first.getEnd(),
										vcAlleles
										);
								vcb.attribute(VCFConstants.END_KEY, first.getEnd());
								vcb.attribute(REPID, getAtt.apply(first, REPID));
								vcb.attribute(RU, getAtt.apply(first, RU));
								vcb.genotypes(genotypes);
								w.add(vcb.make());
							}
						eq.close();
						}
					}
				sorter.cleanup();
				return 0;
				}
			catch(final Throwable err) {
				LOG.error(err);
				return -1;
				}
			finally
				{
			
				}
			}
	
	public static void main(final String[] args)
		{
		new ExpansionHunterMerge().instanceMainWithExit(args);
		}

	}
