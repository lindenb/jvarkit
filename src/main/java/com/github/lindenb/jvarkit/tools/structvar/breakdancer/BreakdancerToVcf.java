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
package com.github.lindenb.jvarkit.tools.structvar.breakdancer;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## input

Input is the tabular output of BreakDancer.

Genotypes are all set to 0/1.

## Example:

```
$ wget -O - -q "https://raw.githubusercontent.com/genome/breakdancer/master/test-data/expected_output" | java -jar dist/breakdancer2vcf.jar -R src/test/resources/human_b37.dict 

##fileformat=VCFv4.2
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimate allele Frequency">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=CHROM2,Number=1,Type=String,Description="Chromosome 2">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=ORIENT1,Number=1,Type=String,Description="Orientation 1">
##INFO=<ID=ORIENT2,Number=1,Type=String,Description="Orientation 2">
##INFO=<ID=POS2,Number=1,Type=String,Description="Position 2">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variation type">
##breakdancer.command=bdfast -o21 inv_del_bam_config
##breakdancer.version=1.4.1-unstable-10-fdfe9f2-dirty (commit fdfe9f2-dirty)
##breakdancer2vcf.meta=compilation:20200511120615 githash:c830e0b htsjdk:2.21.3 date:20200511120941 cmd:-R src/test/resources/human_b37.dict
##contig=<ID=1,length=249250621>
##contig=<ID=2,length=243199373>
##contig=<ID=3,length=198022430>
(...)
##contig=<ID=22,length=51304566>
##contig=<ID=X,length=155270560>
##contig=<ID=Y,length=59373566>
##contig=<ID=MT,length=16569>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H_IJ-NA19238-NA19238-extlibs	H_IJ-NA19240-NA19240-extlibs
21	29185056	.	N	<INS>	99	.	DP=2;END=29185377;SVLEN=226;SVTYPE=INS	GT:DP:GQ	0/1:1:99	0/1:1:99
21	29185462	.	N	<DEL>	99	.	DP=21;END=29186122;SVLEN=-545;SVTYPE=DEL	GT:AF:DP:GQ	0/1:176.58:21:99	0/1:167.89:.:99
21	34807694	.	N	<INS>	99	.	DP=3;END=34808852;SVLEN=304;SVTYPE=INS	GT:DP:GQ	0/1:1:99	0/1:2:99
21	34808937	.	N	<INV>	99	.	DP=2;END=34809799;SVLEN=-737;SVTYPE=INV	GT:AF:DP:GQ	0/1:847.39:.:99	0/1:878.83:2:99
```



END_DOC

 */
@Program(name="breakdancer2vcf",
description="Convert output of breakdancer to VCF",
keywords= {"cnv","sv","breakdancer","vcf"},
creationDate="20200511",
modificationDate="20200511"
)
public class BreakdancerToVcf extends Launcher {
	private static final Logger LOG = Logger.build(BreakdancerToVcf.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-R","-reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path dictRefFile =  null;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate= new WritingVariantsDelegate();

	@Override
	public int doWork(final List<String> args) {
		
		VariantContextWriter out = null;
		try {
			final SAMSequenceDictionary dict = (this.dictRefFile==null?null: SequenceDictionaryUtils.extractRequired(this.dictRefFile));
			
			final Set<VCFHeaderLine> metadata = new HashSet<>();
			
			VCFStandardHeaderLines.addStandardFormatLines(metadata, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.GENOTYPE_QUALITY_KEY,
					VCFConstants.GENOTYPE_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(metadata, true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.END_KEY
					);
			metadata.add(new VCFInfoHeaderLine("SVLEN",1, VCFHeaderLineType.Integer, "SV length"));
			metadata.add(new VCFInfoHeaderLine(VCFConstants.SVTYPE,1, VCFHeaderLineType.String, "Structural variation type"));
			metadata.add(new VCFInfoHeaderLine("ORIENT1",1, VCFHeaderLineType.String, "Orientation 1"));
			metadata.add(new VCFInfoHeaderLine("ORIENT2",1, VCFHeaderLineType.String, "Orientation 2"));
			metadata.add(new VCFInfoHeaderLine("CHROM2",1, VCFHeaderLineType.String, "Chromosome 2"));
			metadata.add(new VCFInfoHeaderLine("POS2",1, VCFHeaderLineType.String, "Position 2"));
						
			final Map<String, String> bam2sample =  new HashMap<>();
			
			
			String header_tokens[]=null;
			try(final BufferedReader br = super.openBufferedReader(oneFileOrNull(args))) {
				String line;
				boolean in_library_statistics = false;
				while((line=br.readLine())!=null) {
					if(!line.startsWith("#")) {
						LOG.error("Bad header line (no hash): "+line);
						return -1;
						}
					else if(line.startsWith("#Software:")) {
						metadata.add(new VCFHeaderLine("breakdancer.version",line.substring(10).trim()));
						}
					else if(line.startsWith("#Command:")) {
						metadata.add(new VCFHeaderLine("breakdancer.command",line.substring(9).trim()));
						}
					else if(line.startsWith("#Library Statistics:")) {
						in_library_statistics=true;
						}
					else if(line.startsWith("#Chr1")) {
						final String stdHeader[] = CharSplitter.TAB.split("#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib");
						header_tokens= CharSplitter.TAB.split(line);
						for(int i=0;i< stdHeader.length ;i++) {
							if(header_tokens.length<=i || !header_tokens[i].equals(stdHeader[i])) {
								LOG.error("Bad header. Expected column $"+(i+1)+"="+stdHeader[i]);
								return -1;
								}
							}
						for(int i=stdHeader.length ;i< header_tokens.length; i++) {
							if(!bam2sample.containsKey(header_tokens[i])) {
								LOG.error("Bam not defined in header ="+header_tokens[i]);
								return -1;
								}
							}
						break;
						}
					else if(in_library_statistics) {
						final String tokens[]= CharSplitter.TAB.split(line);
						for(int i=1;i< tokens.length;i++) {
							final String t=tokens[i];
							if(!t.startsWith("library:")) continue;
							final String sample = t.substring(8);
							if(bam2sample.values().stream().anyMatch(V->V.equals(sample))) {
								LOG.error("Cannot handler Duplicate sample: "+sample);
								return -1;
								}
							String bam = tokens[0].substring(1)/* remove hash */;
							bam2sample.put(bam, sample);
							/* in header , it's the bam name Without the path (!!!) */
							final int slash = bam.lastIndexOf('/');
							if(slash!=-1) {
								bam = bam.substring(slash+1);
								bam2sample.put(bam, sample);
							}
						}
					}
				}
			if(bam2sample.isEmpty()) {
				LOG.error("No sample found in header");
				return -1;
			}
			
			if(header_tokens==null) {
				LOG.error("No header line found.");
				return -1;
			}
			
			final VCFHeader header = new VCFHeader(metadata,new HashSet<>( bam2sample.values()));
			if(dict!=null) header.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			out =  this.writingVariantsDelegate.dictionary(dict).open(this.outputFile);
			out.writeHeader(header);
			while((line=br.readLine())!=null) {
				if(line.startsWith("#")) {
					LOG.error("Bad body "+line);
					return -1;
					}
				final String tokens[]= CharSplitter.TAB.split(line);
				if(tokens.length < 11) {
					LOG.error("Expected 12 columns in "+line.replaceAll("[\t]","<TAB>")+" but got "+tokens.length+". Skipping");
					continue;
					}
				if(dict!=null && dict.getSequence(tokens[0])==null) throw new JvarkitException.ContigNotFoundInDictionary(tokens[0], dict);
				if(dict!=null && dict.getSequence(tokens[3])==null) throw new JvarkitException.ContigNotFoundInDictionary(tokens[3], dict);
				
				final VariantContextBuilder vcb= new VariantContextBuilder();
				vcb.chr(tokens[0]);
				final int start = Integer.parseInt(tokens[1]);
				final int stop = Integer.parseInt(tokens[4]);
				vcb.start(start);
				vcb.attribute(VCFConstants.SVTYPE,tokens[6]);
				
				final Allele REF_ALLELE = Allele.create("N", true);
				final Allele ALT_ALLELE = Allele.create("<" + tokens[6]+">", false);
				vcb.alleles(Arrays.asList(REF_ALLELE,ALT_ALLELE));
				
				
				if(!tokens[0].equals(tokens[3])) {
					vcb.attribute("CHROM2", tokens[3]);
					vcb.attribute("POS2", stop);
					vcb.stop(start);
					}
				else
					{
					vcb.stop(stop);
					vcb.attribute(VCFConstants.END_KEY,stop);
					vcb.attribute("SVLEN", Integer.parseInt(tokens[7])*-1 /* inversed logic */);
					}
				
				vcb.attribute(VCFConstants.DEPTH_KEY, Integer.parseInt(tokens[9]));
				vcb.log10PError(Integer.parseInt(tokens[8])/-10.0);

				
				final List<Genotype> genotypes = new ArrayList<>(bam2sample.size());
				
				for(String bamStr:  CharSplitter.COLON.splitAsStringList(tokens[10]) ) {
						
					final int pipe = bamStr.indexOf('|');
					if(pipe==-1) {
						LOG.error("Pipe missing in "+tokens[10]);
						return -1;
						}
					final String bam = bamStr.substring(0,pipe);
					final int dp = Integer.parseInt(bamStr.substring(pipe+1));

					final String sampleName = bam2sample.get(bam);
					if(StringUtils.isBlank(sampleName)) {
						LOG.error("Cannot get sample associated to bam:"+bam);
						return -1;
						}
					
					final GenotypeBuilder gbuilder = new GenotypeBuilder(sampleName,Arrays.asList(REF_ALLELE,ALT_ALLELE));
					gbuilder.DP(dp);
					gbuilder.GQ(Integer.parseInt(tokens[8]));
					genotypes.add(gbuilder.make());
					}
				vcb.genotypes(genotypes);
				out.add(vcb.make());
				}
			}
			
			out.close();
			out=null;
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(out);
		}
	}

	
	public static void main(final String[] args) {
		new BreakdancerToVcf().instanceMainWithExit(args);
	}

}
