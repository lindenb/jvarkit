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
package com.github.lindenb.jvarkit.tools.vcfgtf;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
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
	modificationDate="20230704",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class VcfFilterGtf extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.of(VcfFilterGtf.class);
	private enum FeatureType {GENE,TRANSCRIPT,EXON,EXON_BOUDARIES};
	private enum Coding {all,protein_coding};
	@Parameter(names={"-gtf","--gtf"},description="GTF file",required = true)
	private Path gtfPath;
	@Parameter(names={"-t","--type"},description="What should be used in the gtf to keep vcf records",required = true)
	private FeatureType featureType = FeatureType.GENE;
	@Parameter(names={"--biotype"},description="What biotype 'x'  should be used.")
	private Coding codingType = Coding.all;

	@Parameter(names={"--extend","--extends","-x"},description="extends each gtf feature by 'x' bases." + DistanceParser.OPT_DESCRIPTION ,splitter = NoSplitter.class, converter = DistanceParser.StringConverter.class)
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
			final List<Locatable> locatables = new ArrayList<>(100_000);
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
					switch(this.featureType) {
						case TRANSCRIPT: ok = feat.getType().equals("transcript"); break;
						case GENE: ok = feat.getType().equals("gene"); break;
						case EXON_BOUDARIES: //continue
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
					if(this.featureType.equals(FeatureType.EXON_BOUDARIES)) {
						for(int side=0;side<2;++side) {
							final int pos = (side==0?feat.getStart():feat.getEnd());
							locatables.add(new SimpleInterval(
								ctg,
								Math.max(1, pos - this.extend),
								pos + this.extend)
								);
							}
						}
					else
						{
						locatables.add(new SimpleInterval(
								ctg,
								Math.max(1, feat.getStart() - this.extend),
								feat.getEnd() + this.extend)
								);
						}
				
					}
 				}
			Collections.sort(locatables,new ContigDictComparator(dict).createLocatableComparator());
			int x=0;
			while(x< locatables.size()) {
				if(x+1< locatables.size() && locatables.get(x).overlaps(locatables.get(x+1))) {
					final Locatable l0 = locatables.get(x);
					final Locatable l1 = locatables.get(x+1);
					locatables.set(x, new SimpleInterval(
							l0.getContig(),
							Math.min(l0.getStart(), l1.getStart()),
							Math.max(l0.getEnd(), l1.getEnd())
							));
					locatables.remove(x+1);
					}
				else
					{
					x++;
					}
				}
			LOG.info("count(intervals)="+locatables.size());

			final IntervalTreeMap<Boolean> intervalTreeMap = new IntervalTreeMap<>();
			for(Locatable loc:locatables) {
				intervalTreeMap.put(new Interval(loc),Boolean.TRUE);
				}
			locatables.clear();
			
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
