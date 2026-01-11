/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.mantamerger;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.sv.StructuralVariantComparator;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

# Input

input is a list of indexed vcf files or one file with the '.list' suffix containing the path to the vcfs


# Example

```
$ find src -type f -name "manta*.vcf.gz" > jeter.list
$ java -jar dist/mantamerger.jar jeter.list 2> /dev/null

(...)
```
END_DOC

 */

@Program(name="mantamerger",
description="Merge Vcf from Manta VCF.",
keywords= {"sv","manta","vcf"},
creationDate="20190916",
modificationDate="20260111",
jvarkit_amalgamion = true,
menu="VCF Manipulation"
)
public class MantaMerger extends Launcher {
	private static final Logger LOG = Logger.of( MantaMerger.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private StructuralVariantComparator svComparator = new StructuralVariantComparator();
	@Parameter(names={"-c","--contig"},description="limit to this contig")
	private String limitContig = null;
	@Parameter(names={"--no-bnd"},description="discar BND")
	private boolean discard_bnd = false;
	@Parameter(names={"-B","--bed"},description="Keep regions in this bed file.")
	private Path includePath = null;

	/**
	 * Minimap description of a genotype: a sample, a genotype type
	 *
	 */
	private static class MiniGT {
		private final String sample;
		private final GenotypeType genotypeType;
		private final boolean passing_filters;
		MiniGT(final String sample,final GenotypeType genotypeType,final boolean passing_filters) {
			this.sample = sample;
			this.genotypeType = genotypeType;
			this.passing_filters = passing_filters;
			}
		
		public String getSample() {
			return sample;
			}
		
		@Override
		public int hashCode() {
			return getSample().hashCode();
			}
		
		public boolean hasSameSample(final  MiniGT o) {
			return getSample().equals(o.getSample());
			}
		@Override
		public boolean equals(final  Object obj) {
			return this.hasSameSample(MiniGT.class.cast(obj));
			}
		private boolean isA(final GenotypeType gt) {
			return this.genotypeType.equals(gt);
			}
		
		boolean isHet() {
			return isA(GenotypeType.HET);
			}
		boolean isHomVar() {
			return isA(GenotypeType.HOM_VAR);
			}
		boolean isHomRefOrNoCall() {
			return !isHet() && !isHomVar();
			}
		}
	/** A unique SV key. A key is defined by a contit, start,reference, and must match the svComparator */
	private class SVKey implements Locatable {
		private final VariantContext archetype;
		SVKey(final VariantContext archetype) {
			this.archetype = archetype;
			}
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof SVKey)) return false;
			final SVKey other = SVKey.class.cast(obj);
			if(!this.contigsMatch(other)) return false;
			if(this.getStart()!=other.getStart()) return false;
			if(!this.archetype.getReference().equals(other.archetype.getReference())) return false;
			if(this.hashCode()!=other.hashCode()) return false;
			return test(other);
			}
		@Override
		public String getContig() {
			return archetype.getContig();
			}
		@Override
		public int getStart() {
			return archetype.getStart();
			}
		@Override
		public int getEnd() {
			return archetype.getEnd();
			}
		public boolean test(SVKey other) {
			return svComparator.test(this.archetype, other.archetype);
			}
		
		@Override
		public int hashCode() {
			int h=1;
			h= h*31 + getContig().hashCode();
			h= h*31 + Objects.hashCode(getStart());
			if(svComparator.isTestingSvTypes()) {
				h= h*31 + this.archetype.getAttributeAsString(VCFConstants.SVTYPE,"").hashCode();
				}
			return h;
			}
		@Override
		public String toString() {
				return new SimpleInterval(this.archetype).toString();
			}
		}
	
	/**
	 * VCF PAth and sample name
	 */
	private static class VcfInput {
		Path vcfPath;
		String sample;
		}

	private boolean isDecoy(final String s) {
		if(s.equals("hs37d5")) return true;
		return false;
	}
	
@Override
public int doWork(final List<String> args) {
	try {
		final Map<String,VcfInput> sample2inputs = new TreeMap<>();
		SAMSequenceDictionary dict=null;
		
		// load input as lines
		final List<String> lines;
		if(args.size()==1 && args.get(0).endsWith(".list"))
			{
			lines= Files.lines(Paths.get(args.get(0))).
				filter(L->!(StringUtils.isBlank(L) || L.startsWith("#"))).
				collect(Collectors.toList())
				;
			}
		else
			{
			lines = args;
			}
		
		/** loop over each input */
		for(final String line: lines) {
			/** split each line tokens[0] is path to VCF, tokens[1] if any is sample */
			final String tokens[]= CharSplitter.TAB.split(line);
			final VcfInput vcfInput = new VcfInput();
			vcfInput.vcfPath  = Paths.get(tokens[0]);
			IOUtil.assertFileIsReadable(vcfInput.vcfPath);
			final SAMSequenceDictionary dict1= SequenceDictionaryUtils.extractRequired(vcfInput.vcfPath);
			if(dict==null) {
				dict=dict1;
				}
			else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict1))
				{
				throw new JvarkitException.DictionariesAreNotTheSame(dict1, dict);
				}
			/* extract sample name from VCF */
			if(tokens.length<2 || StringUtils.isBlank(tokens[1])) {
				try(VCFReader r= VCFReaderFactory.makeDefault().open(vcfInput.vcfPath, false)) {
					final List<String> snl = r.getHeader().getSampleNamesInOrder();
					if(snl.size()==1) {
						vcfInput.sample = snl.get(0);
						}
					else if(snl.isEmpty())
						{
						vcfInput.sample = IOUtils.getFilenameWithoutCommonSuffixes(vcfInput.vcfPath);
						}
					else
						{
						LOG.warn("more than one sample in "+ vcfInput.vcfPath);
						return -1;
						}
					}
				}
			else 
				{
				vcfInput.sample = tokens[1];
				}
			if(sample2inputs.containsKey(vcfInput.sample)) {
				LOG.error("duplicate sample "+vcfInput.sample);
				return -1;
				}
			sample2inputs.put(vcfInput.sample,vcfInput);
			}
			
		if(sample2inputs.isEmpty()) {
			LOG.error("no input found");
			return -1;
			}
		
		if(!StringUtils.isBlank(this.limitContig) && dict.getSequence(this.limitContig)==null) {
			LOG.error(JvarkitException.ContigNotFoundInDictionary.getMessage(this.limitContig, dict));
			return -1;
			}
		
		// load bed containing regions to exclude
		final Predicate<VariantContext> bedPredicate;
		if(this.includePath!=null) {
			final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(dict);
			final IntervalTreeMap<Boolean> map = new IntervalTreeMap<>();
			try(BedLineReader blr = new BedLineReader(this.includePath)) {
				blr.setContigNameConverter(converter);
				blr.stream().
					filter(L->L!=null).
					map(B->new Interval(B)).
					forEach(R->map.put(R,true));
				}
			bedPredicate = V->map.containsOverlapping(V);
			}
		else {
			bedPredicate = V -> true;	
			}
		
		final Set<VCFHeaderLine> metaData = new HashSet<>();
		metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY,true));
		metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
		metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY,true));
		metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY,true));
		metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
		metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY,true));
		metaData.add(new VCFInfoHeaderLine(VCFConstants.SVTYPE, 1, VCFHeaderLineType.String,"Variation type"));

		final VCFInfoHeaderLine infoSvLen = new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"Variation length");
		metaData.add(infoSvLen);

		final VCFInfoHeaderLine infoNSamples = new VCFInfoHeaderLine("NSAMPLES", 1, VCFHeaderLineType.Integer,"Number of samples");
		metaData.add(infoNSamples);
		final VCFInfoHeaderLine infoSamples = new VCFInfoHeaderLine("SAMPLES", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,"amples");
		metaData.add(infoSamples);

		final VCFFilterHeaderLine failingFilters = new VCFFilterHeaderLine("FAIL","All genotypes are FILTERed");
		metaData.add(failingFilters);

		
		
		final VCFFilterHeaderLine filterSameNext = new VCFFilterHeaderLine("NEXT","next variant in VCF is the same.");
		metaData.add(filterSameNext);
		final VCFFilterHeaderLine filterSamePrev = new VCFFilterHeaderLine("PREV","next variant in VCF is the same.");
		metaData.add(filterSamePrev);
		
		
		final VCFHeader header=new VCFHeader(metaData,sample2inputs.keySet()); 
		header.setSequenceDictionary(dict);
		JVarkitVersion.getInstance().addMetaData(this, header);
		
		try(VariantContextWriter out = VCFUtils.createVariantContextWriterToPath(this.outputFile)) {
			out.writeHeader(header);
			//loop over each contig in the dictionary
			for(final SAMSequenceRecord ssr: dict.getSequences()) {
				if(!StringUtils.isBlank(this.limitContig)) {
					if(!ssr.getSequenceName().equals(this.limitContig)) continue;
					}
				
				if(isDecoy(ssr.getSequenceName())) continue;
				LOG.info("contig "+ssr.getSequenceName());

				final List<SVKey> all_uniq_sv_keys = new ArrayList<>(100_000);
				final IntervalTreeMap<Set<SVKey>> intervalTree = new IntervalTreeMap<>();
				
				
				/* scan each input, build a list of intervals */
				for(final VcfInput vcfinput: sample2inputs.values()) {
					/* read all variants in that vcf */
					try(VCFReader vcfFileReader = VCFReaderFactory.makeDefault().open(vcfinput.vcfPath,true)) {
						vcfFileReader.query(ssr).
							stream().
							filter(V->discard_bnd==false || !V.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")).
							filter(bedPredicate).
							map(V->new VariantContextBuilder(V).
									unfiltered().
									noID().
									noGenotypes().
									rmAttribute("EVENT").
									rmAttribute("HOMSEQ").
									rmAttribute("HOMLEN").
									rmAttribute("SVINSSEQ").
									rmAttribute("SVINSLEN").
									rmAttribute("MATEID").
									rmAttribute("LEFT_SVINSSEQ").
									rmAttribute("RIGHT_SVINSSEQ").
									rmAttribute("BND_DEPTH").
									rmAttribute("MATE_BND_DEPTH").
									rmAttribute("JUNCTION_QUAL").
									rmAttribute("CIGAR").
									make()).
							forEach(V->{
								final SVKey key1=new SVKey(V);								
								if(!key1.test(key1)) throw new RuntimeException("compare to self failed ! "+V);
								final Interval r = new Interval(
										key1.getContig(),
										Math.max(1, key1.getStart()- (this.svComparator.getBndDistance() +1)),
										key1.getEnd() + (this.svComparator.getBndDistance() +1)
										);
								
								if(intervalTree.getOverlapping(r).stream().flatMap(coll->coll.stream()).anyMatch(K->K.test(key1))) {
									return;
									}
								Set<SVKey> keys = intervalTree.get(r);
								if(keys==null)  {
									keys = new HashSet<>();
									intervalTree.put(r, keys);
									}
								keys.add(key1);
								all_uniq_sv_keys.add(key1);
							});
						
						}
					}
				if(all_uniq_sv_keys.isEmpty()) continue;
				
	
				// build a map svgkey/genotype
				final Map<SVKey,Set<MiniGT>> variants2genotypes = new HashMap<>(all_uniq_sv_keys.size());
				for(SVKey svkey: all_uniq_sv_keys) {
					variants2genotypes.put(svkey, new HashSet<>());
					}

				// loop over each VCF and scan all genotypes
				for(final VcfInput vcfinput: sample2inputs.values()) {
					try(VCFReader vcfFileReader = VCFReaderFactory.makeDefault().open(vcfinput.vcfPath,true)) {
						try(final CloseableIterator<VariantContext> iter = vcfFileReader.query(ssr)) {
							while(iter.hasNext()) {
								final VariantContext ctx = iter.next();
								
								if(this.discard_bnd && ctx.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")) continue;
								if(!bedPredicate.test(ctx)) continue;
								final SimpleInterval r = new SimpleInterval(ctx).
										extend(this.svComparator.getBndDistance()+1);
								final boolean passing_filters = !ctx.isFiltered();
								final MiniGT miniGt;
								if(ctx.hasGenotypes()) {
									miniGt=new MiniGT(vcfinput.sample,ctx.getGenotype(0).getType(),passing_filters);
									}
								else
									{
									miniGt=new MiniGT(vcfinput.sample,GenotypeType.HET,passing_filters);
									}
								
								
								intervalTree.getOverlapping(r).stream()
									.flatMap(collection->collection.stream())
									.filter(K->this.svComparator.test(K.archetype, ctx))
									.forEach(K->{
			  							final Set<MiniGT> samples= variants2genotypes.get(K);
			  							samples.add(miniGt);
										});
								}
							}
						}
					}
				final Comparator<VariantContext> sorter =new ContigDictComparator(dict).createLocatableComparator();
				final List<SVKey> orderedKeys = variants2genotypes.keySet().
						stream().
						filter(K->!variants2genotypes.get(K).isEmpty()).// no samples for this key ??!
						sorted((A,B)->sorter.compare(A.archetype, B.archetype)).
						collect(Collectors.toCollection(ArrayList::new));
				
				
				
				
				for(int key_index=0;key_index < orderedKeys.size();key_index++) {
					final SVKey key = orderedKeys.get(key_index);
					final Set<MiniGT> miniGts = variants2genotypes.get(key);
					final Allele refAllele = key.archetype.getReference();
					final Object svType = key.archetype.getAttribute(VCFConstants.SVTYPE, ".");
					final Allele altAllele = Allele.create(svType.equals(".")?"<SV>":"<"+svType+">", false);

					final VariantContextBuilder vcb=new VariantContextBuilder();
					vcb.chr(key.getContig());
					vcb.start(key.getStart());
					vcb.stop(key.getEnd());
					vcb.log10PError(key.archetype.getLog10PError());
					vcb.alleles(Arrays.asList(refAllele,altAllele));
					vcb.attribute(VCFConstants.END_KEY, key.getEnd());
					vcb.attribute(VCFConstants.SVTYPE, svType);
					vcb.attribute(infoSvLen.getID(), Math.abs(key.archetype.getLengthOnReference()));
					vcb.attribute(infoNSamples.getID(),miniGts.size());
					vcb.attribute(infoSamples.getID(),miniGts.stream().map(MG->MG.sample).sorted().collect(Collectors.toList()));
					
					int ac = 0;
					final List<Genotype> genotypes = new ArrayList<>(sample2inputs.size());
					for(final String sn:sample2inputs.keySet()) {
						List<Allele> gta;
						final MiniGT miniGt = miniGts.stream().
								filter(MG->MG.sample.equals(sn)).
								findFirst().
								orElse(null);
						if(miniGt==null || miniGt.isHomRefOrNoCall())
							{
							gta =  Arrays.asList(refAllele,refAllele);
							}
						else if(miniGt.isHet()) {
							ac++;
							gta =  Arrays.asList(refAllele,altAllele);
							}
						else if(miniGt.isHomVar()) {
							ac+=2;
							gta =  Arrays.asList(altAllele,altAllele);
							}
						else
							{
							throw new IllegalStateException();
							}
						final GenotypeBuilder gb = new GenotypeBuilder(sn,gta);
						if(miniGt!=null && !miniGt.passing_filters) {
							gb.attribute(VCFConstants.GENOTYPE_FILTER_KEY, "FAIL");
							}
						genotypes.add(gb.make());
						}
					
					if(miniGts.stream().filter(G->!G.isHomRefOrNoCall()).allMatch(G->G.passing_filters==false)) {
						vcb.filter(failingFilters.getID());
						}
					
					vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, ac);
					vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,sample2inputs.size()*2);
					if(ac>0) {
						vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,ac/(float)(sample2inputs.size()*2));
						}
					
					if(key_index>0 && key.test(orderedKeys.get(key_index-1)))
						{
						vcb.filter(filterSamePrev.getID());
						}
					if(key_index+1 < orderedKeys.size()  && key.test(orderedKeys.get(key_index+1)))
						{
						System.err.println("SAME\n"+key.archetype+"\n"+ orderedKeys.get(key_index+1).archetype);
						vcb.filter(filterSameNext.getID());
						}
	
					
					vcb.genotypes(genotypes);
					out.add(vcb.make());
					}
				}
			} //end out
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	}

public static void main(final String[] args) {
	new MantaMerger().instanceMainWithExit(args);
	}

}
