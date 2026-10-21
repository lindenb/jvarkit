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
package com.github.lindenb.jvarkit.tools.vcfgtf;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.LocatableDelegate;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC


## Motivation

filter/annotate a VCF with a GTF. I'm fed up to filter a  GTF, convert the GTF to bed and use the bed to filter the GTF.

END_DOC
*/
@Program(name="vcfgtf",
	description="Filter/Annotate VCF using a GTF file",
	keywords={"vcf","gtf"},
	creationDate="20230703",
	modificationDate="20261021",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class VcfFilterGtf extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.of(VcfFilterGtf.class);
	private enum Coding {all,protein_coding};
	@Parameter(names={"-gtf","--gtf"},description="GTF file",required = true)
	private Path gtfPath;
	@Parameter(names={"--biotype"},description="What biotype 'x'  should be used.")
	private Coding codingType = Coding.all;
	@Parameter(names={"--bed"},description="only use GTF line overlaping that bed file")
	private Path bedPath=null;
	@Parameter(names={"--extend","--extends","-x","-X"},description="extends each gtf feature by 'x' bases. Multiple separated by commas. Empty will be replaced by '0'. " + DistanceParser.OPT_DESCRIPTION)
	private String extendStr = "0";
	@Parameter(names={"--hide-no-overlap","--hide"},description="Don't show variant that have NO intersection with any of the GTF features")
	private boolean hide_no_overlap_flag = false;
	@Parameter(names={"--gene_name"},description="Optional comma delimited string or File containing a list of gene names . Filter on attribute gene_name")
	private String gene_name_path = "";
	@Parameter(names={"--gene_id"},description="Optional comma delimited string or File containing a list of gene IDs . Filter on attribute gene_id. Trailing '.[0-9]+' for version will be removed")
	private String gene_id_path = "";
	@Parameter(names={"--no-gene-id"},description="Do not insert gene ids in the output")
	private boolean hide_gene_ids =false;
	@Parameter(names={"--no-gene-name"},description="Do not insert gene names in the output")
	private boolean hide_gene_names =false;
	@Parameter(names={"--no-biotype"},description="Do not insert gene biotypes in the output")
	private boolean hide_biotypes =false;

	private static class GTFInterval extends LocatableDelegate<GTFLine> {
		protected final int extend;
		GTFInterval(final GTFLine line,int extend) {
			super(line);
			this.extend=extend;
			}
		@Override
		public int getStart() {
			return Math.max(1,getDelegate().getStart()-this.extend);
			}
		@Override
		public int getEnd() {
			return getDelegate().getEnd()+this.extend;
			}
		public String getTag() {
			return "in_"+getDelegate().getType()+ xSuffix(this.extend);
			}
		}
	
	private static class GTFExonBoundaries extends GTFInterval {
		final int side;
		GTFExonBoundaries(GTFLine line,int side,int extend) {
			super(line,extend);
			this.side = side;
			}
		private int getLoc() {
			return this.side==0?getDelegate().getStart():getDelegate().getEnd();
			}
		@Override
		public int getStart() {
			return Math.max(1,getLoc()-this.extend);
			}
		@Override
		public int getEnd() {
			return getLoc()+this.extend;
			}
		@Override
		public String getTag() {
			return "in_exon_boundaries" + xSuffix(this.extend);
			}
		}
	
	/** upstream / downstream gene */
	private static abstract class GTFGeneXXXStream extends GTFInterval {
		protected GTFGeneXXXStream(final GTFLine line,int extend) {
			super(line,extend);
			if(extend<=0) throw new IllegalStateException();
			}
		// each breakend before/after
		protected int x1() { return Math.max(1,getDelegate().getStart()-this.extend);}
		protected int x2() { return Math.max(1,getDelegate().getStart()-1);}
		protected int x3() { return getDelegate().getEnd()+1;}
		protected int x4() { return getDelegate().getEnd()+this.extend;}
		}
	
	/** upstream / downstream gene */
	private static class GTFGeneUpstream extends GTFGeneXXXStream {
		protected GTFGeneUpstream(final GTFLine line,int extend) {
			super(line,extend);
			}
		@Override
		public int getStart() { return getDelegate().isPostiveStrand()?x1():x3();}
		@Override
		public int getEnd() { return getDelegate().isPostiveStrand()?x2():x4();}
		@Override
		public String getTag() {
			return "in_gene_upstream" + xSuffix(this.extend);
			}
		}
	
	/** upstream / downstream gene */
	private static class GTFGeneDownstream extends GTFGeneXXXStream {
		protected GTFGeneDownstream(final GTFLine line,int extend) {
			super(line,extend);
			}
		@Override
		public int getStart() { return getDelegate().isPostiveStrand()?x3():x1();}
		@Override
		public int getEnd() { return getDelegate().isPostiveStrand()?x4():x2();}
		@Override
		public String getTag() {
			return "in_gene_downstream" + xSuffix(this.extend);
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private static String xSuffix(int extend) {
		return  (extend==0?"":"_extended_"+extend);
		}
	
	/** remove ensembl version, if any */
	private static String removeVersion(final String s) {
		final int dot = s.lastIndexOf('.');
		if(dot<=0 || dot+1>=s.length()) return s;
		if(!s.substring(dot+1).matches("[0-9]+")) return s;
		return s.substring(0,dot);
		}
	
	/** load list of identifier (genename, geneid ) from a Path or a string */
	private static Set<String> slurpIdsFromPath(final String stringOrFile,boolean remove_version) throws IOException {
		if(StringUtils.isBlank(stringOrFile)) {
			return null;
			}
		Path as_file =null;
		try {
			as_file =Paths.get(stringOrFile);
			if(!Files.exists(as_file)) as_file = null;
			}
		catch(Throwable err) {
			as_file=null;
			}
		if(as_file==null) {
			return Arrays.stream(stringOrFile.split("[ ,]")).
				map(T->T.trim()).
				map(T->(remove_version?removeVersion(T):T)).
				filter(S->!StringUtils.isBlank(S)).
				collect(Collectors.toSet())
				;
			}
		else
			{
			return Files.readAllLines(as_file).
				stream().
				filter(S->!S.startsWith("#")).
				map(T->T.trim()).
				map(T->(remove_version?removeVersion(T):T)).
				filter(S->!StringUtils.isBlank(S)).
				collect(Collectors.toSet())
				;
			}
		}
	
	private static void pushInterval(final IntervalTreeMap<List<GTFInterval>> treemap,final GTFInterval r) {
		final Interval rgn =new Interval(r);
		List<GTFInterval> L= treemap.get(rgn);
		if(L==null) {
			L=new ArrayList<>();
			treemap.put(rgn, L);
			}
		L.add(r);
		}
	
	private static Set<String> getBioTypes(GTFLine feat) {
		return Arrays.asList(
				feat.getAttribute("biotype"),
				feat.getAttribute("gene_biotype"),
				feat.getAttribute("transcript_biotype")
				)
			.stream()
			.filter(S->!StringUtils.isBlank(S) && !S.equals("."))
			.collect(Collectors.toSet());
		}
	
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator iterin, final VariantContextWriter out) {
		try {
			final Set<Integer> all_extends= CharSplitter.COMMA.stream(this.extendStr)
				.filter(S->!StringUtils.isBlank(S))
				.map(S->new DistanceParser().applyAsInt(S))
				.collect(Collectors.toCollection(TreeSet::new))
				;
			
			if(all_extends.isEmpty()) {
				LOG.warn("extends set is empty, adding '0' as default");
				all_extends.add(0);
				}
			if(all_extends.stream().anyMatch(I->I.intValue()<0)) {
				LOG.error("extends contains a negative value "+this.extendStr);
				return -1;
				}
			
			
			final Set<String> gene_name_set = slurpIdsFromPath(this.gene_name_path,false);
			final Set<String> gene_id_set = slurpIdsFromPath(this.gene_id_path,true);
			final Set<String> allowed_features = new HashSet<>(Arrays.asList(
					"CDS" , "exon"  , "gene" , "start_codon" , "stop_codon" , "transcript" , "UTR",
					"CTCF_binding_site" , "EMAR" , "enhancer" , "open_chromatin_region" , "promoter"
					));
			final VCFHeader header0 = iterin.getHeader();
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header0);
			final UnaryOperator<String> ctgConverter = ContigNameConverter.fromOneDictionary(dict);
			final IntervalTreeMap<List<GTFInterval>> intervalTreeMap = new IntervalTreeMap<>();
			
			
			/* load BED file if needed */
			final IntervalTreeMap<Boolean> bedIntervalMap;
 			if(this.bedPath!=null) {
				try(BedLineReader blr=new BedLineReader(this.bedPath)) {
					blr.setContigNameConverter(ctgConverter);
					bedIntervalMap= blr.toIntervalTreeMap(B->true);
					}
				}
 			else
 				{
 				bedIntervalMap = null;
 				}
			/* load GTF file if needed */
			final GTFCodec codec = new GTFCodec();
			try(BufferedReader br = IOUtils.openPathForBufferedReading(this.gtfPath)) {
				String line;
				while((line=br.readLine())!=null) {
					if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
					final String[] tokens = CharSplitter.TAB.split(line);
					final String ctg = ctgConverter.apply(tokens[0]);
					if(StringUtils.isBlank(ctg)) continue;
					tokens[0]=ctg;
					
					final GTFLine feat = codec.decode(tokens);
					if(feat==null) continue;
					// not in bed
					if(bedIntervalMap!=null && !bedIntervalMap.containsOverlapping(feat)) {
						continue;
						}
					//not a known feature
					if(!allowed_features.contains(feat.getType())) {
						continue;
						}
					
					//check feature type
					boolean ok = false;
					switch(this.codingType) {
						case all: ok=true; break;
						case protein_coding:
							{
							ok= getBioTypes(feat).contains("protein_coding");
							break;
							}
						}
					if(!ok) continue;
					
					/* check gene name */
					if(gene_name_set!=null) {
						final String s= feat.getGeneName();
						if(StringUtils.isBlank(s)) continue;
						if(!gene_name_set.contains(s)) continue;
						}
					/* check gene id */
					if(gene_id_set!=null) {
						final String s= removeVersion(feat.getGeneId());
						if(StringUtils.isBlank(s)) continue;
						if(!gene_id_set.contains(s)) continue;
						}
					
					for(Integer extend: all_extends) {
						pushInterval(intervalTreeMap, new GTFInterval(feat,extend));
						}
					// special features for exon
					if(feat.isExon()) {
						for(int x=0;x<2;++x) {
							for(Integer extend: all_extends) {
								pushInterval(intervalTreeMap, new GTFExonBoundaries(feat,x,extend));
								}
							}
						}
					// special features for genes
					else if(feat.isGene()) {
						for(Integer extend: all_extends) {
							if(extend.intValue()==0) continue;
							GTFGeneXXXStream gxs = new GTFGeneDownstream(feat, extend);
							if(gxs.getStart()< gxs.getEnd()) {
								pushInterval(intervalTreeMap,gxs);
								}
							gxs = new GTFGeneUpstream(feat, extend);
							if(gxs.getStart()< gxs.getEnd()) {
								pushInterval(intervalTreeMap,gxs);
								}
							}
							
						}
	 				}
				}
		
			final VCFHeader header = new VCFHeader(header0);
			JVarkitVersion.getInstance().addMetaData(this, header);
			for(final String s:allowed_features) {
				for(Integer extend: all_extends) {
					final VCFInfoHeaderLine h=new VCFInfoHeaderLine(
							"in_" + s +  xSuffix(extend),
							0,VCFHeaderLineType.Flag,
							"Variant overlap feature "+s+(extend==0?"":" extended by "+extend+" bases")+" in "+this.gtfPath);
					header.addMetaDataLine(h);
					}
				}
			
			for(Integer extend: all_extends) {
				VCFInfoHeaderLine h=new VCFInfoHeaderLine(
						"in_exon_boundaries"+xSuffix(extend),
						0,VCFHeaderLineType.Flag,
						"Variant overlap exon boundaries "+(extend==0?"":" extended by "+extend+" bases")+" in "+this.gtfPath);
				header.addMetaDataLine(h);
				
				if(extend.intValue()>0) {
					h=new VCFInfoHeaderLine(
							"in_gene_downstream"+xSuffix(extend),
							0,VCFHeaderLineType.Flag,
							"Variant downstream gene "+(extend==0?"":" extended by "+extend+" bases")+" in "+this.gtfPath);
					header.addMetaDataLine(h);
					
					h=new VCFInfoHeaderLine(
							"in_gene_upstream"+xSuffix(extend),
							0,VCFHeaderLineType.Flag,
							"Variant upstream gene "+(extend==0?"":" extended by "+extend+" bases")+" in "+this.gtfPath);
					header.addMetaDataLine(h);
					}
				}
			
			final VCFInfoHeaderLine infoGeneNames=new VCFInfoHeaderLine(
					"gene_names",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					"Gene Names from "+this.gtfPath);
			if(!this.hide_gene_names)  {
				header.addMetaDataLine(infoGeneNames);
				}
			
			final VCFInfoHeaderLine infoGeneIds =new VCFInfoHeaderLine(
					"gene_ids",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					"Gene ids from "+this.gtfPath);
			if(!this.hide_gene_ids)  {
				header.addMetaDataLine(infoGeneIds);
				}

			final VCFInfoHeaderLine infoGeneBioTypes =new VCFInfoHeaderLine(
					"biotypes",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
					"Gene/transcripts biotypes from "+this.gtfPath);
			if(!this.hide_biotypes) {
				header.addMetaDataLine(infoGeneBioTypes);
				}
			out.writeHeader(header);
			while(iterin.hasNext()) {
				final VariantContext ctx = iterin.next();
				final List<GTFInterval> intervals = intervalTreeMap.getOverlapping(ctx).stream().flatMap(L->L.stream()).collect(Collectors.toList());
				if(intervals.isEmpty()) {
					if(this.hide_no_overlap_flag) continue;
					out.add(ctx);
					continue;
					}
				
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				Set<String> set = intervals.stream().map(R->R.getTag()).collect(Collectors.toSet());
				
				/* add feature flags VCF */
				for(String tag:set) {
					vcb.attribute(tag, Boolean.TRUE);
					}
				
				/* add gene names in VCF */
				if(!this.hide_gene_names) {
					set = intervals.stream().map(R->R.getDelegate().getGeneName())
							.filter(S->!StringUtils.isBlank(S) && !S.equals("."))
							.collect(Collectors.toCollection(TreeSet::new));
					if(!set.isEmpty()) {
						vcb.attribute(infoGeneNames.getID(),new ArrayList<>(set));
						}
					}
				
				/* add gene ids in VCF */
				if(!this.hide_gene_ids) {
					set = intervals.stream().map(R->R.getDelegate().getGeneId())
							.map(S->removeVersion(S))
							.filter(S->!StringUtils.isBlank(S) && !S.equals("."))
							.collect(Collectors.toCollection(TreeSet::new));
					if(!set.isEmpty()) {
						vcb.attribute(infoGeneIds.getID(),new ArrayList<>(set));
						}
					}
				/* add biotypes in VCF */
				if(!this.hide_biotypes) {
					set = intervals.stream()
							.flatMap(R->getBioTypes(R.getDelegate()).stream())
							.collect(Collectors.toCollection(TreeSet::new));
					
					if(!set.isEmpty()) {
						vcb.attribute(infoGeneBioTypes.getID(),new ArrayList<>(set));
						}
					}

				
				out.add(vcb.make());
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
