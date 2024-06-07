package com.github.lindenb.jvarkit.tools.bed2vcf;
import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Motivation

Convert BED file to VCF, finding REF allele at start and 'N' as ALT allele.
The VCF can be used later for genotyping.

input is a set of BED file or one file with a '.list' extension containing the path to the BED files.

## Example



```
$ bcftools query -f '%CHROM\t%POS0\t%END\n' src/test/resources/rotavirus_rf.vcf.gz | java -jar dist/jvarkit.jar bed2vcf -R src/test/resources/rotavirus_rf.fa | head -n 30
##fileformat=VCFv4.2
##contig=<ID=RF01,length=3302>
##contig=<ID=RF02,length=2687>
##contig=<ID=RF03,length=2592>
##contig=<ID=RF04,length=2362>
##contig=<ID=RF05,length=1579>
##contig=<ID=RF06,length=1356>
##contig=<ID=RF07,length=1074>
##contig=<ID=RF08,length=1059>
##contig=<ID=RF09,length=1062>
##contig=<ID=RF10,length=751>
##contig=<ID=RF11,length=666>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
RF01	970	.	A	N	.	.	.
RF02	251	.	A	N	.	.	.
RF02	578	.	G	N	.	.	.
RF02	877	.	T	N	.	.	.
RF02	1726	.	T	N	.	.	.
RF02	1962	.	TACA	N	.	.	.
RF03	1221	.	C	N	.	.	.
RF03	1242	.	C	N	.	.	.
RF03	1688	.	T	N	.	.	.
RF03	1708	.	G	N	.	.	.
RF03	2150	.	T	N	.	.	.
RF03	2201	.	G	N	.	.	.
RF03	2315	.	G	N	.	.	.
RF03	2573	.	A	N	.	.	.
RF04	887	.	A	N	.	.	.
RF04	991	.	T	N	.	.	.
RF04	1241	.	T	N	.	.	.
```

END_DOC
*/
@Program(name="bed2vcf",
description="Convert BED file to VCF, finding REF allele at start and 'N' as ALT allele",
keywords={"bed","vcf"},
creationDate="20240604",
modificationDate="20240604",
jvarkit_amalgamion = true,
menu="VCF Manipulation"
)
public class BedToVcf extends Launcher {
	private static final Logger LOG = Logger.build(BedToVcf.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	protected Path outputFile=null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path faidx=null;
	@ParametersDelegate
	protected WritingVariantsDelegate writingVariantsDelegate= new WritingVariantsDelegate();
	@Parameter(names={"-s","--symbolic"},description="use symbolic REF allele if bed length > 'x'. 'x' <=0 to disable")
	private int symbolic_length=100;
	
	private void scan(
			final BufferedReader br,
			final VariantContextWriter w,
			final ContigNameConverter ctgConverter,
			final ReferenceSequenceFile reference
			) throws IOException {
		final BedLineCodec codec = new BedLineCodec();
		GenomicSequence genomic = null;
		String line;
		while((line=br.readLine())!=null) {
			if(BedLine.isBedHeader(line)) continue;
			final BedLine rec= codec.decode(line);
			if(rec==null) continue;
			final String ctg  = ctgConverter.apply(rec.getContig());
			if(StringUtils.isBlank(ctg)) continue;
			if(genomic==null || !genomic.hasName(ctg)) {
				genomic = new GenomicSequence(reference, ctg);
				}
			final boolean use_sym_allele = this.symbolic_length >0 && rec.getLengthOnReference()>   this.symbolic_length;
			
			final Allele ref_allele =  use_sym_allele ?
					Allele.create("<REF"+rec.getLengthOnReference()+">", true):
					Allele.create(genomic.subSequence(rec.getBedStart(), rec.getBedEnd()).toString())
					;
			if(ref_allele.equals(Allele.REF_N)) {
				LOG.warn("skipping "+line+" because REF allele is 'N'");
				continue;
				}
			final Allele alt_allele = Allele.ALT_N;
			
			final VariantContextBuilder vcb=new VariantContextBuilder(
					getProgramName(),
					ctg,
					rec.getStart(),
					rec.getEnd(),
					Arrays.asList(ref_allele,alt_allele)
					);
			if(use_sym_allele) {
				vcb.attribute(VCFConstants.END_KEY,rec.getEnd());
				}
			w.add(vcb.make());
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final List<Path> beds = IOUtils.unrollPaths(args);
			try(ReferenceSequenceFile refFile =ReferenceSequenceFileFactory.getReferenceSequenceFile(faidx)) {
				final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(refFile);
				final VCFHeader header= new VCFHeader();
				header.setSequenceDictionary(dict);
				header.addMetaDataLine(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
				final ContigNameConverter ctgConverter = ContigNameConverter.fromOneDictionary(dict);
				try(VariantContextWriter w= this.writingVariantsDelegate.dictionary(dict).open(outputFile)) {
					w.writeHeader(header);
					JVarkitVersion.getInstance().addMetaData(this, header);
					if(beds.isEmpty()) {
						try(BufferedReader br = IOUtils.openStdinForBufferedReader()) {
							scan(br,w,ctgConverter,refFile);
							}
						}
					else
						{
						for(Path path:beds) {
							try(BufferedReader br = IOUtils.openPathForBufferedReading(path)) {
								scan(br,w,ctgConverter,refFile);
								}
							}
						}
					}
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}


	public static void main(String[] args) {
		new BedToVcf().instanceMainWithExit(args);
	}

}
