/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfgtf;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC


## Motivation

filter a VCF with a GTF. I'm fed up to filter a  GTF, convert the GTF to bed and use the bed to filter the GTF.

END_DOC
*/
@Program(name="vcffiltergtf",
	description="Filter VCF on GTF",
	keywords={"vcf","gtf"},
	creationDate="20230703",
	modificationDate="20230703",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class VcfFilterGtf extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VcfFilterGtf.class).make();
	private enum FeatureType {GENE,TRANSCRIPT,EXON};
	private enum Coding {all,protein_coding};
	@Parameter(names={"-gtf","--gtf"},description="GTF file",required = true)
	private Path gtfPath;
	@Parameter(names={"-t","--type"},description="What should be used in the gtf to keep vcf records",required = true)
	private FeatureType featureType = FeatureType.GENE;
	@Parameter(names={"--coding"},description="What type of gtf feature should be used.")
	private Coding codingType = Coding.all;

	@Parameter(names={"--extends"},description="extends each gtf feature by 'x' bases." + DistanceParser.OPT_DESCRIPTION ,splitter = NoSplitter.class, converter = DistanceParser.StringConverter.class)
	private int extend=0;
	@Parameter(names={"--inverse"},description="inverse output")
	private boolean inverse = false;
	@Parameter(names={"--gene_name"},description="Optional file containing a list of gene names . Filter on attribute gene_name")
	private Path gene_name_path = null;
	@Parameter(names={"--gene_id"},description="Optional file containing a list of gene IDs . Filter on attribute gene_id")
	private Path gene_id_path = null;

	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		try {
			final Set<String> gene_name_set;
			final Set<String> gene_id_set;
			
			if(this.gene_name_path!=null) {
				gene_name_set = Files.readAllLines(this.gene_name_path).stream().
					filter(S->!S.startsWith("#")).
					map(T->T.trim()).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet());
				}
			else
				{
				gene_name_set = null;
				}
			
			if(this.gene_id_path!=null) {
				gene_id_set = Files.readAllLines(this.gene_id_path).stream().
					filter(S->!S.startsWith("#")).
					map(T->T.trim()).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet());
				}
			else
				{
				gene_id_set = null;
				}
			
			final VCFHeader header = iterin.getHeader();
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
			final UnaryOperator<String> ctgConverter = ContigNameConverter.fromOneDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, header);
			final IntervalTreeMap<Boolean> intervalTreeMap = new IntervalTreeMap<>();
			final GTFCodec codec = new GTFCodec();
			try(BufferedReader br = IOUtils.openPathForBufferedReading(this.gtfPath)) {
				String line;
				while((line=br.readLine())!=null) {
					if(line.startsWith("#")) continue;
					final GTFLine feat = codec.decode(line);
					if(feat==null) continue;
					final String ctg = ctgConverter.apply(feat.getContig());
					if(StringUtils.isBlank(ctg)) continue;
					boolean ok=false;
					if(!ok) continue;
					switch(this.featureType) {
						case TRANSCRIPT: ok = feat.getType().equals("transcript"); break;
						case GENE: ok = feat.getType().equals("gene"); break;
						case EXON: ok = feat.getType().equals("exon"); break;
						default: break;
						}
					if(!ok) continue;
					
					if(gene_name_set!=null) {
						final String s= feat.getAttribute("gene_name");
						if(StringUtils.isBlank(s)) continue;
						if(!gene_name_set.contains(s)) continue;
						}
					if(gene_id_set!=null) {
						final String s= feat.getAttribute("gene_id");
						if(StringUtils.isBlank(s)) continue;
						if(!gene_id_set.contains(s)) continue;
						}
					
					ok = false;
					switch(this.codingType) {
						case all: ok=true; break;
						case protein_coding:
							{
							String s = feat.getAttribute("transcript_biotype");
							if(s!=null && s.equals("protein_coding")) {
								ok=true;
								break;
								}
							s = feat.getAttribute("gene_biotype");
							if(s!=null && s.equals("protein_coding")) {
								ok=true;
								break;
								}
							s = feat.getAttribute("biotype");
							if(s!=null && s.equals("protein_coding")) {
								ok=true;
								break;
								}
							break;
							}
						}
					if(!ok) continue;
					final Interval r = new Interval(
							ctg,
							Math.max(1, feat.getStart() - this.extend),
							feat.getEnd() + this.extend
							);
					intervalTreeMap.put(r, Boolean.TRUE);
					}
 				}
			out.writeHeader(header);
			while(iterin.hasNext()) {
				final VariantContext ctx = iterin.next();
				boolean ok = intervalTreeMap.containsOverlapping(ctx);
				if(this.inverse) ok=!ok;
				if(!ok) continue;
				out.add(ctx);
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	public static void main(final String[] args)
		{
		new VcfFilterGtf().instanceMainWithExit(args);
		}

}
