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
package com.github.lindenb.jvarkit.tools.upstreamorf;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

Reference: 
	https://github.com/rforge/meet/blob/1ce462dd1732743afd22affdc309d15f7765b657/pkg/R/ModelMATCH.R
	https://www.ncbi.nlm.nih.gov/pubmed/12824369
	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4890304/
	https://github.com/gpnewcomb/uORF_database/blob/034e6a426e55f72a56a068502c5f9df7c0686516/include/defs__general.h#L574
END_DOC

*/
@Program(name="vcfscanupstreamorf",
description="Scan BAM for upstream-ORF",
keywords={"vcf","uorf"},
creationDate="2019-02-08"
)
public class VcfScanUpstreamOrf extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfScanUpstreamOrf.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx = null;
	@Parameter(names={"-k","-K","--kg","-kg"},description=KnownGene.OPT_KNOWNGENE_DESC)
	private String knownGeneUri = KnownGene.getDefaultUri();
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private ContigNameConverter refCtgNameConverter =null;
	private GenomicSequence genomicSequence=null;


	private final IntervalTreeMap<Set<UpstreamOrf>> uOrfMap = new IntervalTreeMap<>();

	private static class KozakBase
		{
		final char base;
		final char choice[];
		KozakBase(char base) {
			this.base= base;
			this.choice = AcidNucleics.degenerateToBases(base);
			}
		boolean match(char c) {
			if(c==this.base) return true;
			for(int i=0;i< choice.length;++i) if(c==choice[i]) return true;
			return false;
			}
		}
	
	private class KozakSequence
		extends AbstractCharSequence
		{
		final List<KozakBase> bases;
		KozakSequence(final List<KozakBase> bases) {
			this.bases=new ArrayList<>(bases);
		}
		KozakSequence(final String bases) {
			this.bases=new ArrayList<>(bases.length());
			for(int i=0;i< bases.length();i++) {
				this.bases.add(new KozakBase(bases.charAt(i)));
			}
		}
		@Override
		public int length() {
			return this.bases.size();
			}
		@Override
		public char charAt(int index) {
			return this.bases.get(index).base;
			}
		int find(final CharSequence dna,int startPos) {
			while(startPos>=0 &&  startPos+this.length() <= dna.length()) {
				// find ATG 
				
				char baseA = Character.toUpperCase(dna.charAt(startPos+0));
				if(baseA!='A') { startPos++; continue;}
				char baseT = Character.toUpperCase(dna.charAt(startPos+1));
				if(baseT!='T') { startPos++; continue;}
				char baseG = Character.toUpperCase(dna.charAt(startPos+2));
				if(baseG!='G') { startPos++; continue;}
				
				return startPos;
				}
			
			return -1;
			}
		}
	
	
	private class UpstreamOrf
		extends AbstractCharSequence
		implements Locatable,Comparable<UpstreamOrf>
		{
		final String contig;
		final int chromStart0;
		final int chromEnd0;
		final boolean plusStrand;
		final Set<String> kgIds = new HashSet<>();
		boolean flag_no_kozak_sequence = false;
		final List<String> matches = null;
	
		private class Match
			extends AbstractCharSequence
			implements Locatable
			{
			int pos=0;
			@Override
			public String getContig() {
				return UpstreamOrf.this.getContig();
				}
			@Override
			public int getStart() {
				return pos;
				}
			@Override
			public int getEnd() {
				return 0;
				}
			@Override
			public int length() {
				return 0;
				}
			@Override
			public char charAt(int index) {
				return 0;
				}
			}
		
		
		UpstreamOrf(final Interval r) {
			this.contig = r.getContig();
			this.chromStart0 = r.getStart()-1;
			this.chromEnd0 = r.getEnd();
			this.plusStrand = !r.isNegativeStrand();
			}
		
		@Override
		public String getContig() {
			return this.contig;
			}
		@Override
		public int getStart() {
			return this.chromStart0 + 1;
			}
		@Override
		public int getEnd() {
			return this.chromEnd0;
			}
		@Override
		public final int length() {
			return this.size();
			}
		public int size() {
			return this.chromEnd0 - chromStart0;
			}
		@Override
		public final char charAt(int index) {
			return readBaseAt0(index);
			}
		public char readBaseAt0(int idx0) {
			if(idx0<0 || idx0>=this.size()) throw new IndexOutOfBoundsException(String.valueOf(idx0+"/size="+size()));
			final char c;
			if(this.plusStrand) {
				final int pos0 = this.chromStart0+idx0;
				c = genomicSequence.charAt(pos0);
				}
			else
				{
				final int pos0 = (this.chromEnd0-1)-idx0;
				c = AcidNucleics.complement(genomicSequence.charAt(pos0));
				}
			return Character.toUpperCase(c);
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof UpstreamOrf)) return false;
			final UpstreamOrf o = UpstreamOrf.class.cast(obj);
			if(o.chromStart0!=this.chromStart0) return false;
			if(o.chromEnd0!=this.chromEnd0) return false;
			if(o.plusStrand!=this.plusStrand) return false;
			return o.contig.equals(this.contig);
			}
		
		@Override
		public int hashCode() {
			int i= this.contig.hashCode();
			i = i*31 + Integer.hashCode(this.chromStart0);
			i = i*31 + Integer.hashCode(this.chromEnd0);
			i = (this.plusStrand?i*31:0);
			return i;
			}
		
		@Override
		public int compareTo(final UpstreamOrf o) {
			int i= contig.compareTo(o.contig);
			if(i!=0)return i;
			i= Integer.compare(chromStart0,o.chromStart0);
			if(i!=0)return i;
			return Integer.compare(chromEnd0,o.chromEnd0);
			}
		}
	
	
		
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator iter,
		final VariantContextWriter out
		) {
		try {
			this.indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidx);
			final SAMSequenceDictionary refDict = SequenceDictionaryUtils.extractRequired(this.indexedFastaSequenceFile);
			this.refCtgNameConverter= ContigNameConverter.fromOneDictionary(refDict);
			
			final KozakSequence kozakSequence = new KozakSequence("ATG");
			
			
			LOG.info("Loading "+this.knownGeneUri);
			try(BufferedReader br= IOUtils.openURIForBufferedReading(this.knownGeneUri)) {
				String line;
				final Map<UpstreamOrf,UpstreamOrf> selfself = new HashMap<>();
				final CharSplitter tab=CharSplitter.TAB;
				while((line=br.readLine())!=null)
					{
					if(StringUtils.isBlank(line))continue;
					final String tokens[]=tab.split(line);
					final String refContig = this.refCtgNameConverter.apply(tokens[0]);
					if(StringUtils.isBlank(refContig)) continue;
					final KnownGene kg=new KnownGene(tokens);
					kg.setChrom(refContig);
					final Interval interval;
					if(kg.isPositiveStrand()) {
						//no utr
						if(kg.getStart()>=kg.getCdsStart()) continue;
						interval = new Interval(kg.getContig(),kg.getStart()+1,kg.getCdsStart(),false,kg.getName());
						}
					else  if(kg.isNegativeStrand())
						{
						//no utr
						if(kg.getCdsEnd()>=kg.getEnd()) continue;
						interval = new Interval(kg.getContig(),kg.getCdsEnd()/* no +1 */,kg.getEnd(),true,kg.getName());
						}
					else
						{
						LOG.info("bad strand in "+kg);
						continue;
						}
					if(interval.length()< kozakSequence.length()) continue;
					
					final UpstreamOrf key = new UpstreamOrf(interval);
					if(selfself.containsKey(key))
						{
						selfself.get(key).kgIds.add(kg.getName());
						}
					else
						{
						key.kgIds.add(kg.getName());
						selfself.put(key, key);
						}
					}
				// end of parsing kg, fill 'this.uOrfMap'
				for(final UpstreamOrf key : selfself.keySet()) {
					final Interval rgn=new Interval(key);
					Set<UpstreamOrf> set =  this.uOrfMap.get(rgn);
					if(set==null) {
						set=new HashSet<UpstreamOrf>();
						this.uOrfMap.put(rgn,set);
						}
					set.add(key);
					}
				LOG.info("number of uORF:"+selfself.size());
				}
  
			if(this.uOrfMap.isEmpty()) {
				LOG.error("no gene found in "+this.knownGeneUri);
				return -1;
				}

			/** build vcf header */
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY,true));
						
			final VCFHeader header= iter.getHeader();
			JVarkitVersion.getInstance().addMetaData(this, header);
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			
			out.writeHeader(header);
			
			while(iter.hasNext()) {
				final VariantContext ctx = progress.apply(iter.next());
				if(!ctx.isVariant()) {
					out.add(ctx);
					continue;
					}
				if(!ctx.isSNP()) {
					out.add(ctx);
					continue;
					}
				final String refContig = this.refCtgNameConverter.apply(ctx.getContig());
				if(StringUtils.isBlank(refContig)) {
					out.add(ctx);
					continue;
					}
				
				/* new reference sequence */
				if(this.genomicSequence==null || !this.genomicSequence.getChrom().equals(refContig)) {
					this.genomicSequence = new GenomicSequence(this.indexedFastaSequenceFile, refContig);
					}

				
				final Interval interval = new Interval(refContig,ctx.getStart(),ctx.getEnd()); 
				final Set<UpstreamOrf> uorfs = this.uOrfMap.getOverlapping(interval).
						stream().
						flatMap(C->C.stream()).
						filter(U-> U.flag_no_kozak_sequence == false).
						collect(Collectors.toSet());
				if(uorfs.isEmpty()) {
					out.add(ctx);
					continue;
					}
				
				for(final UpstreamOrf uorf: uorfs) {
					boolean found_kozak=false;	
					for(int i=0;i+kozakSequence.length()<=uorf.length();++i) {
						int j=0;
						while(j< kozakSequence.length()) {
							if(!kozakSequence.bases.get(j).match(uorf.readBaseAt0(i+j))) {
								break;
							}
						++j;
						}
						
					if(j==kozakSequence.length()) {
						found_kozak = true;
						}
 					}
					
					
	
					if(!found_kozak) {
						uorf.flag_no_kozak_sequence=true;
						}
					}
				
				
				}
			
			progress.close();
			
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args) {
		new VcfScanUpstreamOrf().instanceMainWithExit(args);
	}
}