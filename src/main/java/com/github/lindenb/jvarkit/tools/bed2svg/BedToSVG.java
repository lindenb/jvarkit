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
package com.github.lindenb.jvarkit.tools.bed2svg;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.collections.map.HashedMap;
import org.w3c.dom.Element;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.Hyperlink;
import com.github.lindenb.jvarkit.svg.SVGDocument;
import com.github.lindenb.jvarkit.tools.bed2vcf.BedToVcf;
import com.github.lindenb.jvarkit.tools.misc.CytobandToSvg;
import com.github.lindenb.jvarkit.ucsc.Cytoband;
import com.github.lindenb.jvarkit.util.Maps;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
/**
BEGIN_DOC

## Motivation

Convert BED file to SVG

## Example



```

END_DOC
*/
@Program(name="bed2svg",
description="Convert BED file to VCF, finding REF allele at start and 'N' as ALT allele",
keywords={"bed","svg"},
creationDate="20261013",
modificationDate="20261013",
jvarkit_amalgamion = true,
menu="VCF Manipulation"
)
public class BedToSVG extends Launcher {
	private static final Logger LOG = Logger.of(BedToSVG.class);
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"--cytoband"},description=Cytoband.OPT_DESC)
	private Path cytoBandFile=null;
	@Parameter(names={"-R"},description=DICTIONARY_SOURCE,required = true)
	private Path faidx=null;
	@Parameter(names={"--width"},description="window width in pixels")
	private int screen_width=1000;

	private String ssr_regex = "";
	private long min_ssr_length = 0L;
	
	private final List<Locatable> regions_of_interest = new ArrayList<>();
	private SAMSequenceDictionary dictionary;
	
	private class Track {
		final List<BedLine> records= new ArrayList<>();
		String name="";
		double y=0;
		}
	
	private final Map<String,Track> id2track = new LinkedHashMap<>();
	
	private void readBed(BufferedReader r,final String filename) throws IOException {
		try(BedLineReader blr = new BedLineReader(r,filename)) {
			blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(this.dictionary));
			while(blr.hasNext()) {
				final BedLine rec = blr.next();
				if(!regions_of_interest.isEmpty()) {
					if(regions_of_interest.stream().noneMatch(R->R.overlaps(rec))) continue;
					}
				String trackName= (filename==null?"":filename);
				Track track= id2track.get(trackName);
				if(track==null) {
					track = new  Track();
					track.name = trackName;
					this.id2track.put(trackName, track);
					}
				track.records.add(rec);
				}
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final SAMSequenceDictionary original_dict=SequenceDictionaryUtils.extractRequired(this.faidx);
			final Hyperlink hyperlink = Hyperlink.compile(original_dict);
			dictionary = new SAMSequenceDictionary(
						original_dict
							.getSequences()
							.stream()
							.filter(SR->SR.getSequenceLength() >= this.min_ssr_length)
							.filter(SR->StringUtils.isBlank(ssr_regex) || SR.getContig().matches(ssr_regex))
							.collect(Collectors.toList())
						);
			if(this.dictionary.isEmpty()) {
				LOG.info("no chromosome to display");
				}
			final List<Path> bedPaths= IOUtils.unrollPaths(args);
			
			if(bedPaths.isEmpty()) {
				try(BufferedReader r=new BufferedReader(new InputStreamReader(stdin(), "UTF-8"))) {
					readBed(r,null);
					}
				}
			else
				{
				for(Path bedPath : bedPaths) {
					try(BufferedReader r = IOUtils.openPathForBufferedReading(bedPath)) {
						readBed(r,IOUtils.getFilenameWithoutCommonSuffixes(bedPath));
						}
					}
				}
			
			final SVGDocument svg = new SVGDocument();
			
			double y=0;
			Element g= svg.group();
			final long genome_length = this.dictionary.getReferenceLength();
			Map<String,Long> chromosome2offset=new HashMap<>();
			double featureHeight=10.0;
			long offset = 0L;
			for(SAMSequenceRecord ssr:this.dictionary.getSequences()) {
				chromosome2offset.put(ssr.getContig(), offset);
				Element rect = svg.rect(
					(offset/(double)genome_length)*this.screen_width,
					y,
					(ssr.getLengthOnReference()/(double)genome_length)*this.screen_width,
					featureHeight,
					Maps.of("", "")
					);
				svg.setTitle(rect, ssr.getContig());
				g.appendChild(svg.anchor(rect, hyperlink.apply(ssr).orElse(null)));
				offset+=ssr.getLengthOnReference();
				}
			
			if(this.cytoBandFile!=null) {
				try(BufferedReader br=IOUtils.openPathForBufferedReading(this.cytoBandFile)) {
					final Pattern tab=Pattern.compile("[\t]");
					List<Cytoband> cytobands= br.lines()
						.map(L->tab.split(L))
						.map(T->Cytoband.of(T))
						.filter(C->dictionary.getSequence(C.getContig())!=null)
						.sorted(new ContigDictComparator(this.dictionary).createLocatableComparator())
						.collect(Collectors.toList());
					for(Cytoband cytoband: cytobands) {
						Element rect = svg.rect(
								(cytoband.getStart()/(double)genome_length)*this.screen_width,
								y,
								(cytoband.getLengthOnReference()/(double)genome_length)*this.screen_width,
								featureHeight,
								Maps.of("", "")
								);
						svg.setTitle(rect, cytoband.getName());
						g.appendChild(svg.anchor(rect, hyperlink.apply(cytoband).orElse(null)));
						}
					}
			}
			
			svg.saveTo(this.outputFile);
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	
}
