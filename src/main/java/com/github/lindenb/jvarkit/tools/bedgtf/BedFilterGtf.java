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
package com.github.lindenb.jvarkit.tools.bedgtf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineCodec;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.HasName;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.samtools.util.LocatableDelegate;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;

/**
BEGIN_DOC


## Motivation

filter/annotate a BED with a GTF. Very same idea than  `vcfgtf`


## Example

Filter to protein-coding genes only

```
cat in.bed |\
java -jar jvarkit.jar bedgtf \
    --gtf annotations.gtf.gz \
    --biotype protein_coding \
    --hide-no-overlap \
    --extend 0,1kb \
    -o annotated.bed
```

Annotate structural variants with cut-mode

```
cat in.bed |\
java -jar jvarkit.jar bedgtf \
    --gtf gencode.v44.gtf \
    --cut \
    --extend 0,100 \
    --show-transcript-id \
    --show-protein-id \
    -o sv_annotated.vcf
```

Filter to a subset of genes via file

```
# genes.txt contains one gene name per line
cat in.bed |\
java -jar jvarkit.jar bedgtf \
    --gtf ref.gtf \
    --gene-name genes.txt \
    --hide-no-overlap \
    -o filtered.vcf
```



END_DOC
*/
@Program(name="bedgtf",
	description="Filter/Annotate BED using a GTF file",
	keywords={"vcf","gtf"},
	creationDate="20261113",
	modificationDate="20261113",
	jvarkit_amalgamion =  true
	)
public class BedFilterGtf extends Launcher
	{
	private static final Logger LOG = Logger.of(BedFilterGtf.class);
	private enum Coding {all,protein_coding};
	@Parameter(names={"-gtf","--gtf"},description="GTF file",required = true)
	private Path gtfPath;
	@Parameter(names={"-R","--reference"},description=DICTIONARY_SOURCE,required = true)
	private Path dictSource = null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--biotype"},description="What feature biotype 'x'  should be used.")
	private Coding codingType = Coding.all;
	@Parameter(names={"--bed"},description="only use GTF lines overlaping that bed file. " + BedLineReader.OPT_DESC)
	private Path bedPath=null;
	@Parameter(names={"--extend","--extends","-x","-X"},description="extends each gtf feature by 'x' bases. Multiple separated by commas. Empty will be replaced by '0'. " + DistanceParser.OPT_DESCRIPTION)
	private String extendStr = "0";
	@Parameter(names={"--hide-no-overlap","--hide"},description="Don't show variant that have NO intersection with any of the GTF features")
	private boolean hide_no_overlap_flag = false;
	@Parameter(names={"--gene-name"},description="Optional comma delimited string or File with the '.txt' suffix containing a list of gene names . Filter on attribute gene_name")
	private String gene_name_path = "";
	@Parameter(names={"--gene-id"},description="Optional comma delimited string or File with the '.txt' suffix containing a list of gene IDs . Filter on attribute gene_id. Trailing '.[0-9]+' for version will be removed")
	private String gene_id_path = "";
	@Parameter(names={"--no-gene-id"},description="Do not insert gene ids in the output")
	private boolean hide_gene_ids =false;
	@Parameter(names={"--no-gene-name"},description="Do not insert gene names in the output")
	private boolean hide_gene_names =false;
	@Parameter(names={"--no-biotype"},description="Do not insert gene biotypes in the output")
	private boolean hide_biotypes =false;
	@Parameter(names={"--show-transcript-id"},description="Do show transcript IDs in the output")
	private boolean show_transcripts_id = false;
	@Parameter(names={"--show-protein-id"},description="Do show protein IDs in the output")
	private boolean show_proteins_id = false;
	@Parameter(names={"--cut"},description="The GTF feature should not enclose the variant. This is interesting for the SV like the DUP or the INV where  we want that the variant 'cut' a genomic feature. Features where the whole Variant is contained will be ignored.")
	private boolean cut_flag = false;
	@Parameter(names={"--header-type"},description="Input bed header type. 0 : auto detect ; 1 : the first line is always the header; 2: there is no header.")
	private int input_header_type=0;
	
	private static class Context {
		String line;
		BedLine rec;
		List<GTFInterval> intervals;
		}
	
	private static interface ColumnExtractor extends HasName,Function<Context, String> {
		
		}
	
	
	
	private static class GTFInterval extends LocatableDelegate<GTFLine> {
		protected final int extend;
		GTFInterval(final GTFLine line,final int extend) {
			super(line);
			this.extend=extend;
			}
		protected int trim(int p) {
			return Math.max(1, p);
			}
		@Override
		public int getStart() {
			return trim(getDelegate().getStart()-this.extend);
			}
		@Override
		public int getEnd() {
			return getDelegate().getEnd()+this.extend;
			}
		public String getTag() {
			return "in_"+getDelegate().getType()+ xSuffix(this.extend);
			}
		public String getLabel() {
			String s = getDelegate().getGeneName();
			if(isBlank(s)) s= getDelegate().getTranscriptId();
			if(isBlank(s)) s= getDelegate().getType()+"_"+getContig()+"_"+getStart()+"_"+getEnd();
			return s;
			}
		}
	
	private static class GTFExonBoundaries extends GTFInterval {
		final int side;
		GTFExonBoundaries(final GTFLine line,int side,int extend) {
			super(line,extend);
			this.side = side;
			}
		private int getLoc() {
			return this.side==0?getDelegate().getStart():getDelegate().getEnd();
			}
		@Override
		public int getStart() {
			return trim(getLoc()-this.extend);
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
		protected int x1() { return trim(getDelegate().getStart()-this.extend);}
		protected int x2() { return trim(getDelegate().getStart()-1);}
		protected int x3() { return getDelegate().getEnd()+1;}
		protected int x4() { return getDelegate().getEnd()+this.extend;}
		}
	
	/** upstream  gene */
	private static class GTFGeneUpstream extends GTFGeneXXXStream {
		GTFGeneUpstream(final GTFLine line,int extend) {
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
	
	/**  downstream gene */
	private static class GTFGeneDownstream extends GTFGeneXXXStream {
		GTFGeneDownstream(final GTFLine line,int extend) {
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
	
	private static boolean isBlank(final String s)  {
		return StringUtils.isBlank(s) || s.equals(".");
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
		
		if(stringOrFile.endsWith(".txt")) {
			final Path p = Paths.get(stringOrFile);
			IOUtil.assertFileIsReadable(p);
			return Files.readAllLines(p).
					stream().
					filter(S->!S.startsWith("#")).
					map(T->T.trim()).
					map(T->(remove_version?removeVersion(T):T)).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet())
					;
			}
		else
			{
			return Arrays.stream(stringOrFile.split("[ ,]")).
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
	
	private boolean testOverlap(final BedLine ctx, GTFInterval feature) {
		if(!ctx.overlaps(feature)) throw new IllegalStateException();
		if(this.cut_flag) {
			/* example : a large DUP enclose the variant: reject the variant */
			if(CoordMath.encloses(
					feature.getStart()-1, feature.getEnd()+1,
					ctx.getStart(), ctx.getEnd())
					) {
				return false;
				}
			}
		return true;
		}
	
	private int run(final String inputName, final BufferedReader br,final SAMSequenceDictionary dict) {
		try {
			final String BLANK=".";
			final String INTERNAL_DELIM=",";
			String first_line=br.readLine();
			if(first_line==null) {
				LOG.warn("file "+inputName+" is empty.");
				return 0;
				}
			String[] tokens = CharSplitter.TAB.split(first_line);
			if(tokens.length<3) {
				throw new JvarkitException.TokenErrors(3, tokens);
				}
			final List<String> columns=new ArrayList<>();
			if(this.input_header_type==2 || (this.input_header_type==0 && (dict.getSequence(tokens[0])!=null && 
					StringUtils.isInteger(tokens[1]) &&
					StringUtils.isInteger(tokens[2]) &&
					StringUtils.parseInt(tokens[1]).getAsInt() >=0 && 
					StringUtils.parseInt(tokens[2]).getAsInt() >= StringUtils.parseInt(tokens[1]).getAsInt()
					))) {
				columns.add("chrom");
				columns.add("start");
				columns.add("end");
				for(int i=3;i< tokens.length;i++) {
					columns.add("col."+(i+1));
					}
				}
			else
				{
				for(String h:tokens) columns.add(h);
				first_line = null;
				}
			final List<ColumnExtractor> column_extractors = new ArrayList<>();
			
			
			final Set<Integer> all_extends= CharSplitter.COMMA.stream(this.extendStr)
				.map(S->S.trim())
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
			final UnaryOperator<String> ctgConverter = ContigNameConverter.fromOneDictionary(dict);
			final IntervalTreeMap<List<GTFInterval>> intervalTreeMap = new IntervalTreeMap<>();
			
			
			/* load BED file if needed */
			final IntervalTreeMap<Boolean> bedIntervalMap;
 			if(this.bedPath!=null) {
				try(BedLineReader blr= BedLineReader.openPath(this.bedPath)) {
					blr.setContigNameConverter(ctgConverter);
					bedIntervalMap= blr.toBooleanIntervalTreeMap();
					}
				}
 			else
 				{
 				bedIntervalMap = null;
 				}
			/* load GTF file if needed */
			final GTFCodec gtfCodec = new GTFCodec();
			
			try(BufferedReader br2 = IOUtils.openPathForBufferedReading(this.gtfPath)) {
				String line;
				while((line=br2.readLine())!=null) {
					if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
					tokens = CharSplitter.TAB.split(line);
					final String ctg = ctgConverter.apply(tokens[0]);
					if(StringUtils.isBlank(ctg)) continue;
					tokens[0]=ctg;
					
					final GTFLine feat = gtfCodec.decode(tokens);
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
							ok= feat.getBioTypes().contains("protein_coding");
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
		
			final Set<String> found_tags = intervalTreeMap.values().stream()
					.flatMap(L->L.stream())
					.map(F->F.getTag())
					.collect(Collectors.toSet());

			
			for(final String s:allowed_features) {
				for(final Integer extend: all_extends) {
					final ColumnExtractor c=new ColumnExtractor() {
						public String getName() {
							return "in_" + s +  xSuffix(extend);
							}
						public String apply(final Context ctx) {
							final Set<String> set= ctx.intervals.stream()
								.filter(R->R.getTag().equals(this.getName()))
								.map(R->R.getLabel())
								.collect(Collectors.toSet());
							return set.isEmpty()?BLANK:String.join(INTERNAL_DELIM, set);
							}
						};
					if(found_tags.contains(c.getName())) {
						column_extractors.add(c);
						}
					}
				}
			
			

			
			/* add gene name in BED */
			if(!this.hide_gene_names) {
				column_extractors.add(new ColumnExtractor() {
					public String getName() {
						return "gene_name";
						}
					public String apply(final Context ctx) {
						final Set<String> set = ctx.intervals.stream().map(R->R.getDelegate().getGeneName())
								.filter(S->!isBlank(S))
								.collect(Collectors.toCollection(TreeSet::new));
						return set.isEmpty() ? BLANK : String.join(INTERNAL_DELIM, set);	
						}
					});
			}
			
			/* add gene ids in BED */
			if(!this.hide_gene_ids) {
				column_extractors.add(new ColumnExtractor() {
					public String getName() {
						return "gene_id";
						}
					public String apply(final Context ctx) {
						final Set<String> set = ctx.intervals.stream().map(R->R.getDelegate().getGeneId())
								.map(S->removeVersion(S))
								.filter(S->!isBlank(S))
								.collect(Collectors.toCollection(TreeSet::new));
						return set.isEmpty() ? BLANK : String.join(INTERNAL_DELIM, set);	
						}
					});
				}

			

			/* add biotypes in VCF */
			if(!this.hide_biotypes) {
				column_extractors.add(new ColumnExtractor() {
					public String getName() {
						return "biotype";
						}
					public String apply(final Context ctx) {
						final Set<String> set = ctx.intervals.stream().flatMap(R->R.getDelegate().getBioTypes().stream())
								.filter(S->!isBlank(S))
								.collect(Collectors.toCollection(TreeSet::new));
						return set.isEmpty() ? BLANK : String.join(INTERNAL_DELIM, set);	
						}
					});
				}
			
			
			/* add transcript ids in BED */
			if(this.show_transcripts_id) {
				column_extractors.add(new ColumnExtractor() {
					public String getName() {
						return "transcript_id";
						}
					public String apply(final Context ctx) {
						final Set<String> set = ctx.intervals.stream().map(R->R.getDelegate().getTranscriptId())
								.map(S->removeVersion(S))
								.filter(S->!isBlank(S))
								.collect(Collectors.toCollection(TreeSet::new));
						return set.isEmpty() ? BLANK : String.join(INTERNAL_DELIM, set);	
						}
					});
				}
			/* add protein ids in VCF */
			if(this.show_proteins_id) {
				column_extractors.add(new ColumnExtractor() {
					public String getName() {
						return "protein_id";
						}
					public String apply(final Context ctx) {
						final Set<String> set = ctx.intervals.stream().map(R->R.getDelegate().getProteinId())
								.map(S->removeVersion(S))
								.filter(S->!isBlank(S))
								.collect(Collectors.toCollection(TreeSet::new));
						return set.isEmpty() ? BLANK : String.join(INTERNAL_DELIM, set);	
						}
					});
				}
			
			final BedLineCodec codec = new BedLineCodec();
			final Context context = new Context();
			codec.setContigNameConverter(dict);
			
			try(PrintWriter out=super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				for(int i=0;i< columns.size();i++) {
					if(i>0) out.print("\t");
					out.print(columns.get(i));
					}
				for(int i=0;i< column_extractors.size();i++) {
					out.print("\t");
					out.print(column_extractors.get(i).getName());
					}
				out.println();
				for(;;) {
					if(first_line!=null) {
						context.line = first_line;
						first_line=null;
						}
					else
						{
						context.line = br.readLine();
						}
					if(context.line==null) break;
					if(StringUtils.isBlank(context.line)) continue;
					tokens = CharSplitter.TAB.split(context.line);
					if(tokens.length > columns.size()) {
						throw new JvarkitException.TokenErrors(columns.size(), tokens);
						}
					context.rec = codec.decode(context.line);
					if(context.rec==null) {
						LOG.warn("skipping  "+context.line);
						continue;
						}
					context.intervals = intervalTreeMap
							.getOverlapping(context.rec)
							.stream()
							.flatMap(L->L.stream())
							.filter(feat->testOverlap(context.rec,feat))
							.collect(Collectors.toList());
					
					if(context.intervals.isEmpty() && this.hide_no_overlap_flag) {
						continue;
						}
					
					out.print(context.rec.toBed3());
					for(int i=3;i< columns.size();i++) {
						out.print("\t");
						out.print(context.rec.getOrDefault(i, ""));
						}
					for(int i=0;i< column_extractors.size();i++) {
						out.print("\t");
						out.print(column_extractors.get(i).apply(context));
						}
					out.println();					
					}
				out.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	@Override
	public int doWork(final List<String> args) {
		try {
			if(this.input_header_type<0 || this.input_header_type>2) {
				LOG.error("bad value for --header-type (0|1|2) but got "+this.input_header_type);
				return -1;
				}
			final String input = oneAndOnlyOneFile(args);
			final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(this.dictSource);
			try(BufferedReader br= input==null?
					IOUtils.openStreamForBufferedReader(stdin()):
					IOUtils.openURIForBufferedReading(input)) {
			
				return run((input==null?"<STDIN>":input),br,dict);
				}
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	public static void main(final String[] args)
		{
		new BedFilterGtf().instanceMainWithExit(args);
		}

}
