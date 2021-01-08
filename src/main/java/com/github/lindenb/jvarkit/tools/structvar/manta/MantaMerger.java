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
package com.github.lindenb.jvarkit.tools.structvar.manta;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.Decoy;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.sv.StructuralVariantComparator;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
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
 $ find src -name "manta*z" > jeter.list
 $ java -jar dist/mantamerger.jar jeter.list 2> /dev/null
 
 (...)
 
 
 
 ```
 END_DOC

 */

@Program(name="mantamerger",
description="Merge Vcf from Manta VCF.",
keywords= {"sv","burden","manta","vcf"},
creationDate="20190916",
modificationDate="20191106"
)
public class MantaMerger extends Launcher {
	private static final Logger LOG = Logger.build( MantaMerger.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@ParametersDelegate
	private StructuralVariantComparator svComparator = new StructuralVariantComparator();
	@Parameter(names={"-c","--contig"},description="limit to this contig")
	private String limitContig = null;
	@Parameter(names={"--no-bnd"},description="discar BND")
	private boolean discard_bnd = false;
	@Parameter(names={"-x","--exclude"},description="Exclude regions in this bed file")
	private Path excludeBedPath = null;

	
	private class SVKey {
		final VariantContext archetype;
		SVKey(final VariantContext archetype) {
			this.archetype = archetype;
			}
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof SVKey)) return false;
			final SVKey other = SVKey.class.cast(obj);
			if(!this.archetype.contigsMatch(other.archetype)) return false;
			if(this.archetype.getStart()!=other.archetype.getStart()) return false;
			if(!this.archetype.getReference().equals(other.archetype.getReference())) return false;
			if(this.hashCode()!=other.hashCode()) return false;
			return svComparator.test(this.archetype, other.archetype);
			}
		@Override
		public int hashCode() {
			int h=1;
			h= h*31 + archetype.getContig().hashCode();
			h= h*31 + archetype.getStart();
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
	

	private static class VcfInput {
		Path vcfPath;
		String sample;
		int contigCount=0;
		}
	
@Override
public int doWork(final List<String> args) {
	VariantContextWriter out = null;
	try {
		final Map<String,VcfInput> sample2inputs = new TreeMap<>();
		SAMSequenceDictionary dict=null;
		
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
		
			
		for(final String line: lines) {
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
			if(tokens.length<2 || StringUtils.isBlank(tokens[1])) {
				try(VCFReader r= VCFReaderFactory.makeDefault().open(vcfInput.vcfPath, false)) {
					List<String> snl = r.getHeader().getSampleNamesInOrder();
					if(snl.size()==1) {
						vcfInput.sample = snl.get(0);
						}
					else
						{
						vcfInput.sample = vcfInput.vcfPath.toString();
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
		
		
		final Predicate<VariantContext> bedPredicate;
		if(this.excludeBedPath!=null) {
			final BedLineCodec codec = new BedLineCodec();
			final ContigNameConverter converter = ContigNameConverter.fromOneDictionary(dict);
			final IntervalTreeMap<Boolean> map = new IntervalTreeMap<>();
			try(BufferedReader br=IOUtils.openPathForBufferedReading(this.excludeBedPath))
				{
				br.lines().
					filter(L->!BedLine.isBedHeader(L)).
					map(L->codec.decode(L)).
					filter(L->L!=null).
					filter(B->!StringUtils.isBlank(converter.apply(B.getContig()))).
					map(B->new Interval(converter.apply(B.getContig()), B.getStart(),B.getEnd())).
					forEach(R->map.put(R,true));
				}
			bedPredicate = V->!map.containsOverlapping(V);
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
		metaData.add(new VCFInfoHeaderLine(VCFConstants.SVTYPE, 1, VCFHeaderLineType.String,"Variation type"));

		final VCFInfoHeaderLine infoSvLen = new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"Variation length");
		metaData.add(infoSvLen);

		final VCFInfoHeaderLine infoNSamples = new VCFInfoHeaderLine("NSAMPLES", 1, VCFHeaderLineType.Integer,"Number of samples");
		metaData.add(infoNSamples);
		final VCFInfoHeaderLine infoSamples = new VCFInfoHeaderLine("SAMPLES", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,"amples");
		metaData.add(infoSamples);

		
		final VCFFilterHeaderLine filterSameNext = new VCFFilterHeaderLine("NEXT","next variant in VCF is the same.");
		metaData.add(filterSameNext);
		final VCFFilterHeaderLine filterSamePrev = new VCFFilterHeaderLine("PREV","next variant in VCF is the same.");
		metaData.add(filterSamePrev);
		
		
		final VCFHeader header=new VCFHeader(metaData,sample2inputs.keySet()); 
		header.setSequenceDictionary(dict);
		JVarkitVersion.getInstance().addMetaData(this, header);
		
		out = VCFUtils.createVariantContextWriterToPath(this.outputFile);
		out.writeHeader(header);

		
		
		
		final Decoy decoy = Decoy.getDefaultInstance();
		for(final SAMSequenceRecord ssr: dict.getSequences()) {
			if(!StringUtils.isBlank(this.limitContig)) {
				if(!ssr.getSequenceName().equals(this.limitContig)) continue;
				}
			
			LOG.info("contig "+ssr.getSequenceName());
			if(decoy.isDecoy(ssr.getSequenceName())) continue;
			final Map<SVKey,Set<String>> variants2samples = new HashMap<>();
			for(final VcfInput vcfinput: sample2inputs.values()) {
				vcfinput.contigCount=0;//reset count for this contig
				
				try(VCFReader vcfFileReader = VCFReaderFactory.makeDefault().open(vcfinput.vcfPath,true)) {
					vcfFileReader.query(ssr.getSequenceName(), 1, ssr.getSequenceLength()).
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
							if(!svComparator.test(V, V)) throw new RuntimeException("compare to self failed ! "+V);
							variants2samples.put(key1, new HashSet<>());
							vcfinput.contigCount++;
						});
					
					}
				}
			if(variants2samples.isEmpty()) continue;
			
			// build an interval tree for a faster access
			final IntervalTree<SVKey> intervalTree = new IntervalTree<>();
			for(final SVKey key: variants2samples.keySet()) {
				final SimpleInterval r = new SimpleInterval(key.archetype).
						extend(this.svComparator.getBndDistance()+1);
				intervalTree.put(r.getStart(), r.getEnd(), key);
				
				// paranoidcheck interval is ok to find current archetype
				boolean found=false;
				final Iterator<IntervalTree.Node<SVKey>> nodeIter = intervalTree.overlappers(r.getStart(),r.getEnd());
				while( nodeIter.hasNext()) {
					final SVKey key1 = nodeIter.next().getValue();
					if(this.svComparator.test(key1.archetype, key.archetype))  {
						found=true;
						break;
						}
					}
				if(!found) {
					out.close();
					throw new RuntimeException("cannot find self "+key.archetype+" in "+r);
					}
				}

			
			for(final VcfInput vcfinput: sample2inputs.values()) {
				try(VCFReader vcfFileReader = VCFReaderFactory.makeDefault().open(vcfinput.vcfPath,true)) {
					
					final CloseableIterator<VariantContext> iter = vcfFileReader.query(ssr.getSequenceName(), 1, ssr.getSequenceLength());
					while(iter.hasNext()) {
						final VariantContext ctx = iter.next();
						
						if(this.discard_bnd && ctx.getAttributeAsString(VCFConstants.SVTYPE, "").equals("BND")) continue;
						if(!bedPredicate.test(ctx)) continue;
						final SimpleInterval r = new SimpleInterval(ctx).
								extend(this.svComparator.getBndDistance()+1);
						
						
						final Iterator<IntervalTree.Node<SVKey>> nodeIter = intervalTree.overlappers(r.getStart(),r.getEnd());
						while( nodeIter.hasNext()) {
							final SVKey key1 = nodeIter.next().getValue();
							if(!this.svComparator.test(key1.archetype, ctx)) continue;
  							final Set<String> samples= variants2samples.get(key1);
  							samples.add(vcfinput.sample);
							}
						}
					iter.close();
					}
				}
			final Comparator<VariantContext> sorter =new ContigDictComparator(dict).createLocatableComparator();
			final List<SVKey> orderedKeys = variants2samples.keySet().
					stream().
					filter(K->!variants2samples.get(K).isEmpty()).// no samples for this key ??!
					sorted((A,B)->sorter.compare(A.archetype, B.archetype)).
					collect(Collectors.toCollection(ArrayList::new));
			
			// remove duplicates
			int i=0;
			while(i +1 < orderedKeys.size())
				{
				final SVKey key1 = orderedKeys.get(i);
				final SVKey key2 = orderedKeys.get(i+1);
				if(svComparator.test(key1.archetype,key2.archetype) &&
					variants2samples.get(key1).equals(variants2samples.get(key2)) // same samples
					) {
					orderedKeys.remove(i+1);
					}
				else
					{
					i++;
					}
				}
			
			
			for(int key_index=0;key_index < orderedKeys.size();key_index++) {
				final SVKey key = orderedKeys.get(key_index);
				final Set<String> samples = variants2samples.get(key);
				final Allele refAllele = key.archetype.getReference();
				final Allele altAllele = Allele.create("<SV>", false);
				final Object svType = key.archetype.getAttribute(VCFConstants.SVTYPE, ".");
				final VariantContextBuilder vcb=new VariantContextBuilder();
				vcb.chr(key.archetype.getContig());
				vcb.start(key.archetype.getStart());
				vcb.stop(key.archetype.getEnd());
				vcb.log10PError(key.archetype.getLog10PError());
				vcb.alleles(Arrays.asList(refAllele,altAllele));
				vcb.attribute(VCFConstants.END_KEY, key.archetype.getEnd());
				vcb.attribute(VCFConstants.SVTYPE, svType);
				vcb.attribute(infoSvLen.getID(), (svType.equals("DEL")?-1:1)*key.archetype.getLengthOnReference());
				vcb.attribute(infoNSamples.getID(),samples.size());
				vcb.attribute(infoSamples.getID(),samples.stream().sorted().collect(Collectors.toList()));
				
				int ac = 0;
				final List<Genotype> genotypes = new ArrayList<>(sample2inputs.size());
				for(final String sn:sample2inputs.keySet()) {
					List<Allele> gta;
					if(samples.contains(sn))
						{
						ac++;
						gta =  Arrays.asList(refAllele,altAllele);
						}
					else
						{
						gta =  Arrays.asList(refAllele,refAllele);
						}
					genotypes.add(new GenotypeBuilder(sn,gta).make());
					}
				
				vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, ac);
				vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,sample2inputs.size()*2);
				if(ac>0) {
					vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,ac/(float)sample2inputs.size()*2);
					}
				
				if(key_index>0 && svComparator.test(key.archetype,  orderedKeys.get(key_index-1).archetype))
					{
					vcb.filter(filterSamePrev.getID());
					}
				if(key_index+1 < orderedKeys.size()  && svComparator.test(key.archetype,  orderedKeys.get(key_index+1).archetype))
					{
					System.err.println("SAME\n"+key.archetype+"\n"+ orderedKeys.get(key_index+1).archetype);
					vcb.filter(filterSameNext.getID());
					}

				
				vcb.genotypes(genotypes);
				out.add(vcb.make());
				}
			}
		out.close();
		out=null;
		return 0;
		}
	catch(final Throwable err) {
		LOG.error(err);
		return -1;
		}
	finally
		{
		CloserUtil.close(out);
		}
	}

public static void main(final String[] args) {
	new MantaMerger().instanceMainWithExit(args);
	}

}
