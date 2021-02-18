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
package com.github.lindenb.jvarkit.tools.biostar;

import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFReader;


/**
BEGIN_DOC

## History

major bug detected after https://github.com/lindenb/jvarkit/issues/178 . Some reads with INDEL might have been ignored.

##Example

```
# look at the original bam at position 14:79838836
$ samtools tview -d T -p 14:79838836 src/test/resources/HG02260.transloc.chr9.14.bam | cut -c 1-20
     79838841  79838
NNNNNNNNNNNNNNNNNNNN
CATAGGAAAACTAAAGGCAA
cataggaaaactaaaggcaa
cataggaaaaCTAAAGGCAA
cataggaaaactaaaggcaa
CATAGGAAAACTAAAGGCAA
cataggaaaactaaaggcaa
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
           taaaggcaa

# look at the VCF , we want to change 14:79838836 with 'A' with probability = 0.5
$ cat jeter.vcf
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency among genotypes, for each ALT allele, in the same order as listed">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
14	79838836	.	N	A	.	.	AF=0.5

# apply biostar404363
$  java -jar dist/biostar404363.jar -o jeter.bam -p jeter.vcf src/test/resources/HG02260.transloc.chr9.14.bam

# index the sequence
$ samtools index jeter.bam

# view the result (half of the bases were replaced because AF=0.5 in the VCF)
$ samtools tview -d T -p 14:79838836 jeter.bam | cut -c 1-20
     79838841  79838
NNNNNNNNNNNNNNNNNNNN
MATAGGAAAACTAAAGGCAA
cataggaaaactaaaggcaa
aataggaaaaCTAAAGGCAA
cataggaaaactaaaggcaa
CATAGGAAAACTAAAGGCAA
aataggaaaactaaaggcaa
AATAGGAAAACTAAAGGCAA
AATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
CATAGGAAAACTAAAGGCAA
AATAGGAAAACTAAAGGCAA
           taaaggcaa
```

END_DOC
*/
@Program(name="biostar404363",
	keywords={"sam","bam","variant"},
	description="introduce artificial mutation SNV in bam",
	biostars= {404363,416897},
	creationDate="20191023",
	modificationDate="20191024"
	)
public class Biostar404363 extends Launcher {
	
	private static final Logger LOG = Logger.build(Biostar404363.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--position","--vcf"},description="VCF File containing the positions to change. if INFO/AF(allele frequency) field is present, variant is inserted if rand()<= AF. ",required=true)
	private Path vcfPath=null;
	@Parameter(names= {"-R","--reference"},description=CRAM_INDEXED_REFENCE)
	private Path faidx;
	@Parameter(names= {"--disable-nm"},description="Disable NM change. By default the value of NM is increasing to each change.")
	private boolean keep_original_nm = false;
	@Parameter(names= {"--ignore-clip"},description="Ignore clipped section of the reads.")
	private boolean ignore_clip = false;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	private static class PosToChange extends SimplePosition {
		byte base = 'N';
		double proba=1.0;
		PosToChange(final String ctg,int pos) {
			super(ctg,pos);
		}
		@Override
		public String toString() {
			return super.toString()+" "+(char)base+" "+proba;
			}
	}

	@Override
	public int doWork(final List<String> args) {
		SamReader in=null;
		SAMFileWriter out=null;
		final IntervalTreeMap<PosToChange> intervalTreeMap = new IntervalTreeMap<>();
		try {
			
			final SamReaderFactory srf = super.createSamReaderFactory();
			if(this.faidx!=null) {
				srf.referenceSequence(this.faidx);
				writingBamArgs.setReferencePath(this.faidx);
				}
			final String input = oneFileOrNull(args);
			if(input==null) {
				in = srf.open(SamInputResource.of(stdin()));
				}
			else
				{
				in = srf.open(SamInputResource.of(input));
				}
			final SAMFileHeader header = in.getFileHeader();
			final SAMProgramRecord prg =  header.createProgramRecord();
			prg.setCommandLine(this.getProgramCommandLine());
			prg.setProgramName(this.getProgramName());
			prg.setProgramVersion(this.getGitHash());
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final ContigNameConverter contigNameConverter=ContigNameConverter.fromOneDictionary(dict);
			
			try(VCFReader vcfReader= VCFReaderFactory.makeDefault().open(this.vcfPath,false)) {
				try(CloseableIterator<VariantContext> iter=vcfReader.iterator()) {
					while(iter.hasNext()) {
						final VariantContext ctx=iter.next();
						final List<Allele> alts = ctx.getAlternateAlleles();
						if(alts.size()!=1) {
							LOG.error("Expected one ALT allele in "+ctx);
							return -1;
							}
						final Allele alt = alts.get(0);
						if(alt.isSymbolic() || alt.length()!=1 || !AcidNucleics.isATGCN(alt)) {
							LOG.error("Bad ALT allele in "+ctx);
							return -1;
							}
						final String theContig= contigNameConverter.apply(ctx.getContig());
						if(StringUtils.isBlank(theContig)) throw new JvarkitException.ContigNotFoundInDictionary(ctx.getContig(), dict);

						final PosToChange pos2change = new PosToChange(theContig, ctx.getStart());
						pos2change.base=(byte)Character.toUpperCase(alt.getBases()[0]);
						if(ctx.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
							pos2change.proba =  ctx.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY,1.0);
							}	
						final Interval key=new Interval(pos2change);
						if(intervalTreeMap.containsKey(key)) {
							LOG.error("Duplicate key "+key+ " for "+ctx);
							return -1;
							}
						intervalTreeMap.put(key,pos2change);
					}
				}
			}
			
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			out = this.writingBamArgs.openSamWriter(this.outputFile, header, true);
			final SAMRecordIterator iter=in.iterator();
			while(iter.hasNext()) {
				final SAMRecord rec= iter.next();
				final byte bases[]=rec.getReadBases();

				if(rec.getReadUnmappedFlag() ||
					bases==SAMRecord.NULL_SEQUENCE) {
					out.addAlignment(rec);
					continue;
					}
				
				final List<PosToChange> changes = intervalTreeMap.getOverlapping(new SimpleInterval(rec.getContig(), rec.getUnclippedStart(), rec.getUnclippedEnd())).
						stream().
						collect(Collectors.toList());
				
				if(changes.isEmpty()) {
					out.addAlignment(rec);
					continue;
					}
				
				
				int NM=0;
				if(rec.hasAttribute("NM")) NM=rec.getIntegerAttribute("NM");
				int readpos=0;
				int ref=rec.getUnclippedStart();
				boolean changed=false;
				final Cigar cigar = rec.getCigar();
				
				for(final CigarElement ce:cigar) {
					final CigarOperator op =ce.getOperator();
					if( (op.equals(CigarOperator.S) && !this.ignore_clip) || 
						(op.consumesReadBases() &&  op.consumesReferenceBases() )) {
						for(int i=0;i< ce.getLength();i++) {
							final int the_pos = ref+i;
							final byte the_base = bases[readpos+i];
							final PosToChange pos2change = changes.stream().
									filter(P->P.getPosition()==the_pos).
									filter(P->P.base!=the_base).
									filter(P->Math.random()<=P.proba).
									findFirst().
									orElse(null);
							if(pos2change==null) continue;
							bases[readpos+i]=pos2change.base;
							if(op.isAlignment()) NM++;
							changed=true;
							}
						}
					if(op.equals(CigarOperator.S) || op.consumesReadBases()) {
						readpos+=ce.getLength();
						}
					if(op.isClipping() || op.consumesReferenceBases()) {
						ref+=ce.getLength();
						}
					}
				
				if(changed) {
					rec.setReadBases(bases);
					if(!keep_original_nm) rec.setAttribute("NM", NM);
					rec.setAttribute("PG", prg.getId());
					}
				out.addAlignment(rec);
				}
			in.close();
			in=null;
			out.close();
			out=null;
			
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(in);
			CloserUtil.close(out);
		}
	}
	
	public static void main(final String[] args) {
		new Biostar404363().instanceMainWithExit(args);
	}

}
