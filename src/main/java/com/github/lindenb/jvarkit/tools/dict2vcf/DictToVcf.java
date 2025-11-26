/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.dict2vcf;

import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.VcfHeaderExtractor;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## Motivation

convert Dict to VCF when you want to generate an empty vcf file

## Example

```
java -jar dist/jvarkit.jar dict2vcf src/test/resources/human_b37.dict  --samples S1,S3
```


END_DOC
 */
@Program(name="dict2vcf",
	description="convert a SAM dictionary from vcf,sam,bam,dict, etc.. to vcf.",
	keywords={"dict","bed","sam","bam","vcf"},
	creationDate="20251126",
	modificationDate="20251126",
	jvarkit_amalgamion =  true
	)
public class DictToVcf extends Launcher {
	private static Logger LOG=Logger.of(DictToVcf.class);

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"--info","--infos"},description="Add those predefined INFO fields")
	private String info_fiels = String.join(",", VCFConstants.END_KEY,VCFConstants.DBSNP_KEY,VCFConstants.DEPTH_KEY,VCFConstants.STRAND_BIAS_KEY,VCFConstants.ALLELE_FREQUENCY_KEY,VCFConstants.ALLELE_COUNT_KEY,VCFConstants.ALLELE_NUMBER_KEY,VCFConstants.MAPPING_QUALITY_ZERO_KEY,VCFConstants.RMS_MAPPING_QUALITY_KEY,VCFConstants.SOMATIC_KEY);
	@Parameter(names={"--format","--formats"},description="Add those predefined FORMAT fields")
	private String formats_fiels = String.join(",", VCFConstants.GENOTYPE_KEY,VCFConstants.GENOTYPE_QUALITY_KEY,VCFConstants.DEPTH_KEY,VCFConstants.GENOTYPE_PL_KEY,VCFConstants.GENOTYPE_ALLELE_DEPTHS,VCFConstants.GENOTYPE_FILTER_KEY,VCFConstants.PHASE_SET_KEY,VCFConstants.PHASE_QUALITY_KEY);
	@Parameter(names={"--extra-lines"},description="File containing extra header lines")
	private Path extraLines=null;
	@Parameter(names={"--samples"},description="Add those samples")
	private List<String> samples= new ArrayList<>();
	@Parameter(names={"--samples-file"},description="Read samples from file. If the file has a BAM/CRAM/BCF/VCF extension, samples are extracted from their metadata.")
	private Path samples_file = null;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();
	
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final String input = super.oneFileOrNull(args);
			final SAMSequenceDictionary dict;
			if(input==null  || input.equals("-")) {
				final VCFHeader h = VcfHeaderExtractor.decode(stdin());
				dict = SequenceDictionaryUtils.extractRequired(h);
				}
			else
				{
				dict = SequenceDictionaryUtils.extractRequired(Paths.get(input));
				}
			final Set<VCFHeaderLine> metadata = new HashSet<>();
			
			VCFStandardHeaderLines.addStandardInfoLines(metadata, true,
					Arrays.stream(this.info_fiels.split("[,; ]+"))
						.filter(S->!StringUtils.isBlank(S))
						.collect(Collectors.toSet())
					);
			VCFStandardHeaderLines.addStandardFormatLines(metadata, true,
					Arrays.stream(this.formats_fiels.split("[,; ]+"))
						.filter(S->!StringUtils.isBlank(S))
						.collect(Collectors.toSet())
					);
			Set<String> samples = new TreeSet<>(this.samples.stream().
					flatMap(S->Arrays.stream(S.split("[,; ]+"))).
						filter(S->!StringUtils.isBlank(S)).
						collect(Collectors.toSet())
					);
			if(samples_file!=null) {
				final String fname=samples_file.getFileName().toString();
				
				if(fname.endsWith(FileExtensions.BAM) || fname.endsWith(FileExtensions.CRAM) ) {
					SamReaderFactory.make().getFileHeader(samples_file).getReadGroups().stream()
						.map(RG->RG.getSample())
						.filter(S->!StringUtils.isBlank(S))
						.forEach(S->samples.add(S));
					}
				else if(FileExtensions.VCF_LIST.stream().anyMatch(SUFF->fname.endsWith(SUFF))) {
					samples.addAll(VcfHeaderExtractor.decode(samples_file).getSampleNamesInOrder());
					}
				else
					{
					Files.readAllLines(this.samples_file).stream().
						filter(S->!StringUtils.isBlank(S))
						.forEach(S->samples.add(S));
					}
				}
			
			final VCFHeader h=new VCFHeader(metadata,samples);
			JVarkitVersion.getInstance().addMetaData(this, h);
			if(extraLines!=null) {
				final VCFHeaderVersion version = h.getVCFHeaderVersion();
				for(String str: Files.readAllLines(extraLines, Charset.defaultCharset())) {
					if(StringUtils.isBlank(str)) continue;
					if ( str.startsWith(VCFConstants.INFO_HEADER_START) ) {
	                    h.addMetaDataLine(new VCFInfoHeaderLine(str.substring(7), version));
	                } else if ( str.startsWith(VCFConstants.FILTER_HEADER_START) ) {
	                	h.addMetaDataLine(new VCFFilterHeaderLine(str.substring(9), version));
	                } else if ( str.startsWith(VCFConstants.FORMAT_HEADER_START) ) {
	                	h.addMetaDataLine(new VCFFormatHeaderLine(str.substring(9), version));
	                } else
	                {
	                	LOG.warn("Cannot convert line "+str+" to VCF line header. May be not supported");
	                	return -1;
	                }
				}
			}
			
			h.setSequenceDictionary(dict);
		
			
			try(VariantContextWriter w=writingVariants.dictionary(dict).open(this.outputFile)) {
				w.writeHeader(h);
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new DictToVcf().instanceMainWithExit(args);
	}

}
