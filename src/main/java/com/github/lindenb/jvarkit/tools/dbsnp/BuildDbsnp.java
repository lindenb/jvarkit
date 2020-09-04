/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.dbsnp;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

/**
BEGIN_DOC

## Example

```
$ java  -jar dist/depthofcoverage.jar -R src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam 2> /dev/null  | column -t 

#BAM                       Sample  Contig  Length  Count   Depth
src/test/resources/S1.bam  S1      RF01    3302    25037   7.582374318594791
src/test/resources/S1.bam  S1      RF02    2687    20275   7.545589877186453
src/test/resources/S1.bam  S1      RF03    2592    19583   7.55516975308642
src/test/resources/S1.bam  S1      RF04    2362    17898   7.577476714648603
src/test/resources/S1.bam  S1      RF05    1579    11887   7.528182393920202
src/test/resources/S1.bam  S1      RF06    1356    10201   7.522861356932153
src/test/resources/S1.bam  S1      RF07    1074    8115    7.555865921787709
src/test/resources/S1.bam  S1      RF08    1059    7980    7.5354107648725215
src/test/resources/S1.bam  S1      RF09    1062    7980    7.5141242937853105
src/test/resources/S1.bam  S1      RF10    751     5740    7.6431424766977365
src/test/resources/S1.bam  S1      RF11    666     5037    7.563063063063063
src/test/resources/S1.bam  S1      *       18490   139733  7.557220118983234
src/test/resources/S2.bam  S2      RF01    3302    25030   7.580254391278014
src/test/resources/S2.bam  S2      RF02    2687    20272   7.544473390398213
src/test/resources/S2.bam  S2      RF03    2592    19592   7.558641975308642
src/test/resources/S2.bam  S2      RF04    2362    17916   7.585097375105843
src/test/resources/S2.bam  S2      RF05    1579    11892   7.531348955034832
src/test/resources/S2.bam  S2      RF06    1356    10217   7.534660766961652
src/test/resources/S2.bam  S2      RF07    1074    8112    7.553072625698324

```

END_DOC
 */
@Program(name="builddbsnp",
	description="Build a DBSNP file from different sources for GATK",
	keywords={"vcf","dbsnp"},
	creationDate="20200904",
	modificationDate="20200904",
	generate_doc=false
	)
public class BuildDbsnp extends Launcher {
	private static Logger LOG=Logger.build(BuildDbsnp.class).make();

	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-c","--chromosome"},description="limit to this chromosome")
	private String limitChrom = null;

	
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	
	private static class Variant implements Comparable<Variant> {
		int pos;
		List<Allele> alleles;
		String id;
		@Override
		public int compareTo(final Variant o) {
			int i = Integer.compare(this.pos, o.pos);
			if(i!=0) return i;
			i = alleles.get(0).compareTo(o.alleles.get(0));
			return i;
			}
	}
	
	private static class VCFSource extends AbstractIterator<Variant> implements CloseableIterator<Variant> {
		String name;
		Path filePath;
		VCFReader reader;
		String currentContig;
		CloseableIterator<VariantContext> iter;
		final List<Variant> stack  = new ArrayList<>();
		VCFSource(final String name,final Path path) {
			this.name = name;
			this.filePath = path;
			
			}
		void reset(final SAMSequenceRecord ssr) {
			close();
			this.reader = VCFReaderFactory.makeDefault().open(this.filePath, true);
			this.currentContig = ssr.getContig();
			final VCFHeader header =  this.reader.getHeader();
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict!=null) {
				final String vcfCtg = ContigNameConverter.fromOneDictionary(dict).apply(ssr.getContig());
				if(StringUtils.isBlank(vcfCtg)) {
					iter = this.reader.query(ssr.getContig()/* returns an empty iterator*/,ssr.getStart(),ssr.getEnd());
					}
				else
					{
					iter = this.reader.query(vcfCtg,ssr.getStart(),ssr.getEnd());
					}
				}
			else
				{
				final List<String> contigs = new ArrayList<>();
				contigs.add(ssr.getContig());
				if(ssr.getSequenceName().toLowerCase().startsWith("chr")) {
					contigs.add(ssr.getContig().substring(3));
					}
				else
					{
					contigs.add("chr"+ssr.getContig());
					contigs.add("CHR"+ssr.getContig());
					}
				if(ssr.getSequenceName().equals("chrM") || ssr.getSequenceName().equals("chrMT") ||
					ssr.getSequenceName().equals("M") || ssr.getSequenceName().equals("MT")) {
					contigs.add("chrMT");
					contigs.add("chrM");
					contigs.add("MT");
					contigs.add("M");
					}
				for(final String contig : contigs) {
					iter = this.reader.query(contig,ssr.getStart(),ssr.getEnd());
					if(iter.hasNext()) break;
					iter.close();
					}
				}
			}
		
		@Override
		protected Variant advance() {
			if(!stack.isEmpty()) {
				return stack.remove(0);
				}
			while(iter.hasNext()) {
				final VariantContext ctx = iter.next();
				final Variant variant = new Variant();
				variant.id = ctx.hasID()?ctx.getID():this.name+"_"+this.currentContig+"_"+ctx.getStart();
				variant.alleles = ctx.getAlleles();
				variant.pos = ctx.getStart();
				stack.add(variant);
				if(stack.size()>1) {
					final Variant first = stack.get(0);
					if(first.pos!=variant.pos) break;
					}
				}
			if(!stack.isEmpty()) {
				if(stack.size()>1) Collections.sort(this.stack);
				return stack.remove(0);
				}
			return null;
			}
		
		@Override
		public void close()  {
			if(this.iter!=null) iter.close();
			try { this.reader.close();} catch(Throwable err) {} 
			}
	}
	
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.faidx);
			if(!StringUtils.isBlank(this.limitChrom)) {
				if(dict.getSequence(this.limitChrom)==null) {
					throw new JvarkitException.ContigNotFoundInDictionary(this.limitChrom, dict);
				}
			}
			
			final List<VCFSource> sources = new ArrayList<>(args.size());
			
			for(int i=0;i< args.size();i++) {
				final String s = args.get(i);
				final int colon = s.indexOf(':');
				if(colon<=0) {
					LOG.error("colon missing in "+s+" expected vcf-name:vcf-path.");
					return -1;
					}
				final String vcfName = s.substring(0,colon);
				if(sources.stream().anyMatch(S->S.name.equals(vcfName)))  {
					LOG.error("duplicate source name: \""+vcfName+"\".");
					return -1;
					}
				final Path path = Paths.get(s.substring(colon+1));
				IOUtil.assertFileIsReadable(path);
				sources.add(new VCFSource(vcfName, path));
				}
			final Pattern rsIdPattern = Pattern.compile("[rR][sS][0-9]+");
			try(VariantContextWriter w=writingVariantsDelegate.dictionary(dict).open(this.outputFile)) {
				final VCFHeader header = new VCFHeader();
				header.setSequenceDictionary(dict);
				JVarkitVersion.getInstance().addMetaData(this, header);
				w.writeHeader(header);
				for(final SAMSequenceRecord ssr:dict.getSequences()) {
					if(!StringUtils.isBlank(this.limitChrom) && !ssr.getContig().equals(this.limitChrom)) continue;
					sources.stream().forEach(SRC->SRC.reset(ssr));
					final List<CloseableIterator<Variant>> iterators  = sources.stream().map(SRC->SRC).collect(Collectors.toList());
					final MergingIterator<Variant> merger = new MergingIterator<>((A, B)->A.compareTo(B),iterators);
					final EqualRangeIterator<Variant> equal_range = new EqualRangeIterator<>(merger, (A, B)->A.compareTo(B));
					
					while(equal_range.hasNext()) {
						final List<Variant> variants = equal_range.next();
						final Variant variant = variants.get(0);
						final Set<Allele> altSet = new HashSet<>();
						variants.stream().forEach(V->altSet.addAll(V.alleles.subList(1, V.alleles.size())));
						altSet.remove(Allele.SPAN_DEL);
						final List<Allele> alleles = new ArrayList<>(1+altSet.size());
						alleles.add(variant.alleles.get(0));
						alleles.addAll(altSet);
						
						
						final  VariantContextBuilder vcb = new VariantContextBuilder(null, ssr.getContig(), variant.pos,variant.pos+variant.alleles.get(0).length()-1, alleles);
						String id = variants.stream().filter(F->rsIdPattern.matcher(F.id).matches()).map(F->F.id).findFirst().orElse(null);
						if(StringUtils.isBlank(id)) id = variants.stream().map(F->F.id).findFirst().orElse(null);
						vcb.id(variant.id);
						w.add(vcb.make());
						}
				
					equal_range.close();
					merger.close();
					iterators.stream().forEach(R->R.close());
					}
				}
			sources.stream().forEach(SRC->SRC.close());
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		}
	}
	
	public static void main(final String[] args) {
		new BuildDbsnp().instanceMainWithExit(args);
	}

}
