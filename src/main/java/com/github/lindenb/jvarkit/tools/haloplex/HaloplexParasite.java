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
package com.github.lindenb.jvarkit.tools.haloplex;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**

BEGIN_DOC

### Examples

```

echo "input.bam" > all.list
gunzip -c input.vcf.gz |
  java -jar dist/haloplexparasite.jar -B all.list
rm all.list


```


END_DOC
*/
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC

### Examples

```

echo "input.bam" > all.list
gunzip -c input.vcf.gz |
  java -jar dist/haloplexparasite.jar -B all.list
rm all.list



```




END_DOC
 */
@Program(name="haloplexparasite",
	description="for @SolenaLS : remove artctifacts from haloplex that gives indels in GATK hapcaller ",
	keywords={"vcf","haloplex"}
	)
public class HaloplexParasite extends Launcher {
	private static final Logger LOG = Logger.build(HaloplexParasite.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-B","--bams"},description="A list of path to indexed BAM files")
	private File bamList = null;

	@Parameter(names={"-m","--clipsize"},description="Min. Soft Clipping size")
	private int minClipSize = 10 ;

	@Parameter(names={"-t","--treshold"},description="treshold")
	private double tresholdFraction = 0.0005 ;

	private static class Mutation {
		final String contig;
		final int start;
		final int end;
		final List<Allele> alleles;
		Mutation(final VariantContext ctx) {
			this.contig = ctx.getContig();
			this.start = ctx.getStart();
			this.end = ctx.getEnd();
			this.alleles = ctx.getAlleles();
		}
	}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator in, VariantContextWriter out) {
	SamReader samReader = null;
	final List<Mutation> mutations = new ArrayList<>();
	try {
		final Set<File> bamFiles =  Files.lines(this.bamList.toPath()).
			filter(T-> !(T.isEmpty() || T.startsWith("#"))).
			map(T->new File(T)).collect(Collectors.toSet());
		
		final VCFHeader header = new VCFHeader();
		header.setSequenceDictionary(in.getHeader().getSequenceDictionary());
		final VCFFilterHeaderLine filter = new VCFFilterHeaderLine("HALOPLEXPARASITE", "fails Parasite Haloplex Sequence");
		final VCFInfoHeaderLine infoWord = new VCFInfoHeaderLine(filter.getID(),1,VCFHeaderLineType.String, "Parasite Haloplex Sequence (Word|Count|Fraction)");
		
		super.addMetaData(header);
		out.writeHeader(header);
		header.addMetaDataLine(filter);
		header.addMetaDataLine(infoWord);
		
		while(in.hasNext()) {
			final VariantContext ctx = in.next();
			final VariantContextBuilder vcb = new VariantContextBuilder(
					inputName,
					ctx.getContig(),
					ctx.getStart(),
					ctx.getEnd(),
					ctx.getAlleles()
					);
			
			if(!(ctx.isIndel() ||ctx.isMixed())) {
				out.add(vcb.make());
				continue;
			}
			if( !vcb.getAlleles().stream().
					filter(A->! (A.isSymbolic() || A.isNoCall() || A.length() < this.minClipSize )).
					findAny().isPresent()) {
				out.add(vcb.make());
				continue;
			}
			
		final Mutation mut = new Mutation(ctx);
		mutations.add(mut);
		}
		
		final Counter<String> words = new Counter<>();
		for(final File bamFile : bamFiles) {
			LOG.info("Opening "+bamFile);
			samReader = createSamReaderFactory().open(bamFile);
			for(final Mutation mut: mutations) {
			//words seen in that mutation : don't use a Counter
			final Set<String> mutWords = new HashSet<>();
			/* loop over reads overlapping that mutation */
			final SAMRecordIterator sri = samReader.queryOverlapping(mut.contig,mut.start,mut.end);
			while(sri.hasNext()) {
				final SAMRecord rec = sri.next();
				if(rec.getReadUnmappedFlag() ) continue;
				if(rec.isSecondaryOrSupplementary()) continue;
				if(rec.getDuplicateReadFlag()) continue;
				if(rec.getReadFailsVendorQualityCheckFlag()) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar.numCigarElements()==1) continue;
				final byte bases[]=rec.getReadBases();
				int refpos = rec.getUnclippedStart();
				int readpos = 0;
				/* scan cigar overlapping that mutation */
				for(final CigarElement ce: cigar) {
					final CigarOperator op =ce.getOperator();
					final int ref_end = refpos + (op.consumesReferenceBases() || op.isClipping() ? ce.getLength() : 0 );
					final int read_end = readpos + (op.consumesReadBases()? ce.getLength() : 0 );
					/* check clip is large enough */
					if(op.equals(CigarOperator.S) && ce.getLength()>=this.minClipSize) {
						/* check clip overlap mutation */
						if(!(ref_end < mut.start ||  refpos > mut.end ) ) {
							/* break read of soft clip into words */
							for(int x=0;x + this.minClipSize <= ce.getLength();++x ) {
								final String substr = new String(bases, readpos+ x, this.minClipSize);
								if(!substr.contains("N")) {
									final String revcomp = AcidNucleics.reverseComplement(substr);
									mutWords.add(substr);
									if(!revcomp.equals(substr)) mutWords.add(revcomp);
									}
							}
						}
					}
					refpos = ref_end;
					readpos = read_end;
					}
				}	
			sri.close();
			
			for(final String w: mutWords) {
			words.incr(w);
			}
			
			}//end of for(mutations)
			samReader.close();
			samReader=null;
		}
		
		LOG.info("mutations:"+mutations.size());
		
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);
		for(final Mutation mut: mutations) {
			progress.watch(mut.contig, mut.start);
			final VariantContextBuilder vcb = new VariantContextBuilder(
					inputName,
					mut.contig,
					mut.start,
					mut.end,
					mut.alleles
					);
			
			String worstString = null;
			Double worstFraction=null;
			final double totalWords = words.getTotal(); 
			for(final Allele a: mut.alleles) {
				if(a.isSymbolic() || a.isNoCall() || a.length()< this.minClipSize) continue;
				for(int x = 0;x+this.minClipSize <= a.length();++x) {
					final String substr = new String(a.getBases(),x, this.minClipSize);
					final long count = words.count(substr);
					final double fraction=count/totalWords;
					if(worstFraction==null || worstFraction < fraction) {
						worstString=substr+"|"+count+"|"+fraction;
						worstFraction = fraction;
					}
				}
			}
			if(worstString!=null) {
				vcb.attribute(infoWord.getID(), worstString);
			}
			
			if(worstFraction!=null && worstFraction.doubleValue()>=this.tresholdFraction) {
				vcb.filter(filter.getID());
				
			}
			
			out.add(vcb.make());
		}
		progress.finish();
		
		return RETURN_OK;
	} catch (Exception e) {
		LOG.error(e);
		return -1;
	} finally {
		CloserUtil.close(samReader);
		}
	}
	@Override
	public int doWork(List<String> args) {
		if(this.bamList==null) {
			LOG.error("BAM list undefined");
			return -1;
		}
		return doVcfToVcf(args,outputFile);
		}
	
	public static void main(final String[] args) {
		new HaloplexParasite().instanceMainWithExit(args);
	}
}
