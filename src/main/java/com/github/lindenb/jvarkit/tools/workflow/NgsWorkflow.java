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
package com.github.lindenb.jvarkit.tools.workflow;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonPrimitive;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
/*
BEGIN_DOC

## Examples

```
{
"name":"TSV",
"output.directory":"/path/toTSV",
"prefix":"20170531.TSV.",
"use.lumpyexpress":true,
"samples":[
{"name":"S1","bam":"S1.reliable.realign.bam"},
{"name":"S2","bam":"S2.reliable.realign.bam"}
]
}

```
END_DOC

*/
@Program(
		name="ngsworkflow",
		description="ngs workflow",
		keywords={"ngs","workflow","pipeline","bam","vcf"}
		)
public class NgsWorkflow extends Launcher
	{
	private static final Logger LOG =Logger.build(NgsWorkflow.class).make();
	
	private enum RefSplitType {WHOLE_GENOME,WHOLE_CONTIG,INTERVAL};
	
	@Parameter(names={"-A","--attributes"},description="Dump available attributes and exit")
	private boolean dumpAttributes = false;

	@Parameter(names={"--reannotatevcf"},description="This project is just reannotating vcfs")
	private boolean reannotatevcfs = false;

	
	private abstract class RefSplit
		{
		abstract RefSplitType getType();
		abstract String getToken();
		}
	
	private class NoSplit extends RefSplit
		{
		@Override RefSplitType getType() { return RefSplitType.WHOLE_GENOME;}
		@Override String getToken() { return "";}
		}
	
	private class ContigSplit extends RefSplit
		{
		final String contig;
		ContigSplit(final String contig) {this.contig=contig;}
		public String getContig() {return contig;}
		@Override RefSplitType getType() { return RefSplitType.WHOLE_CONTIG;}
		@Override String getToken() { return "."+getContig();}
		}

	
	private class IntervalSplit extends RefSplit
		{
		final Interval interval;
		IntervalSplit(final Interval interval) {this.interval=interval;
			if(interval.getStart()>=interval.getEnd()) throw new IllegalArgumentException(interval.toString());
			}
		IntervalSplit(String c,int S,int E) {this(new Interval(c,S,E));}
		IntervalSplit(int c,int S,int E) {this(String.valueOf(c),S,E);}

		public Interval getInterval() {return interval;}
		@Override RefSplitType getType() { return RefSplitType.INTERVAL;}
		@Override String getToken() { return "."+getInterval().getContig()+"_"+getInterval().getStart()+"_"+getInterval().getEnd();}
		}
		
	private static interface PropertyKey
		{
		public String getKey();
		public String getDescription();
		public String getDefaultValue();
		}
	
	private static class PropertyKeyBuilder
		{
		private final PropertyKeyImpl prop;
		PropertyKeyBuilder(final String key)
			{
			if(NgsWorkflow.name2propertyKey.containsKey(key))
				{
				throw new IllegalStateException(key);
				}
			this.prop = new PropertyKeyImpl(key);
			}
		PropertyKeyBuilder description(final String s) { this.prop.description=s;return this;}
		PropertyKeyBuilder def(final String s) { this.prop.def=s;return this;}
		PropertyKeyBuilder def(final boolean b) {return def(String.valueOf(b));}
		PropertyKeyBuilder def(final int b) {return def(String.valueOf(b));}
		
		PropertyKey build() {
			NgsWorkflow.name2propertyKey.put(this.prop.key,this.prop);
			return this.prop;
			}
		private  static class PropertyKeyImpl implements PropertyKey
			{
			final String key;
			String description="No Description available";
			String def;
			PropertyKeyImpl(final String key) {
				this.key=key;
				
			}
			
			public @Override String getKey() { return this.key;}
			public @Override String getDescription() { return this.description;}
			public @Override String getDefaultValue() { return this.def;}
			
			@Override
			public int hashCode() {
				return key.hashCode();
				}
			@Override
			public boolean equals(final Object obj) {
				if(obj==this) return true;
				if(obj==null || !(obj instanceof PropertyKeyImpl)) return false;
				return this.key.equals(PropertyKeyImpl.class.cast(obj).key);
				}
			@Override
			public String toString() {
				return this.key.toString();
				}
			}
		
		}
	
	private static final Map<String,PropertyKey> name2propertyKey=new HashMap<>();
	private static PropertyKeyBuilder key(final String key) { return new PropertyKeyBuilder(key);}
	
	
	private static final PropertyKey PROP_PROJECT_OUTPUT_DIRECTORY = key("output.directory").
			description("Project output directory").
			build();
	
	private static final PropertyKey PROP_TMP_PREFIX_TOKEN = key("tmp.prefix").
			description("Temporary Files prefix").
			def("__DELETE__").
			build();
	
	private static final PropertyKey PROP_FILES_PREFIX_TOKEN = key("prefix").
		description("Files prefix").
		build();
	
	private static final PropertyKey PROP_CUTADAPT_ADAPTER5 = key("cutadapt.adapter5").
		description("Cutadapt 5' adapter").
		def("CTACACTCTTTCCCTACACGACGCTCTTCCGATCT").
		build();
	
	private static final PropertyKey PROP_CUTADAPT_ADAPTER3 = key("cutadapt.adapter3").
		description("Cutadapt 3' adapter").
		def("GATCGGAAGAGCACACGTCTGAACTCCAGTCAC").
		build();
	
	private static final PropertyKey PROP_BWA_MEM_NTHREADS = key("bwa.mem.nthreads").
			description("Number of threads for bwa mem").
			def("1").
			build();
	
	private static final PropertyKey PROP_SAMTOOLS_SORT_NTHREADS = key("samtools.sort.nthreads").
			description("Number of threads for samtools sort").
			def(1).
			build();

	private static final PropertyKey PROP_PICARD_MARKDUP_REMOVES_DUPLICATES = key("markdups.removes.duplicates").
			description("picard MarkDuplicate will remove the dups").
			def(false).
			build();
	
	private static final PropertyKey PROP_CAPTURE_EXTEND = key("capture.extend").
			description("How to extend the capture bed").
			def(1000).
			build();

	private static final PropertyKey PROP_DEFAULT_COMPRESSION_LEVEL = key("compression.level").
			description("Compression level").
			def(9).
			build();
	
	private static final PropertyKey PROP_DISABLE_RECALIBRATION = key("disable.recalibration").
			description("Disable GATK Recalibration").
			def(false).
			build();
	private static final PropertyKey PROP_DISABLE_REALIGN = key("disable.realign").
			description("Disable GATK RealignAroundIndel").
			def(false).
			build();
	private static final PropertyKey PROP_DISABLE_MARKDUP = key("disable.markdup").
			description("Disable picard MarkDups").
			def(false).
			build();
	private static final PropertyKey PROP_PROJECT_IS_HALOPLEX = key("haloplex").
			description("Project is Haloplex").
			def(false).
			build();
	private static final PropertyKey PROP_DISABLE_CUTADAPT = key("disable.cutadapt").
			description("Disable cutadapt").
			def(false).
			build();
	private static final PropertyKey PROP_GATK_HAPCALLER_NCT = key("gatk.hapcaller.nct").
			description("Hapcaller -nct").
			def("1").
			build();

	private static final PropertyKey PROP_MAPPING_REGION = key("mapping.region").
			description("When mapping with bwa, use bedtools intersect to restrict the output to the specified bed file").
			def("").
			build();
	private static final PropertyKey PROP_USE_LUMPY_EXPRESS = key("use.lumpyexpress").
			description("Use Lumpy Express").
			def("false").
			build();

	
	private RefSplit parseRefSplitFromStr(final String s)
		{
		int colon= s.indexOf(":");
		if(colon!=-1)
			{
			int hyphen = s.lastIndexOf('-', colon+1);
			if(hyphen<0) throw new IllegalArgumentException("Illegal array item "+s);
			final Interval interval = new Interval(
					s.substring(0,colon),
					Integer.parseInt(s.substring(colon+1,hyphen)),
					Integer.parseInt(s.substring(hyphen+1))
					);
			return new IntervalSplit(interval);
			}
		else
			{
			return new ContigSplit(s);
			}
		}
	
	private List<RefSplit> parseRefSplitList(final JsonElement root) throws IOException
		{
		final List<RefSplit> splits=new ArrayList<>();
		if(root.isJsonObject())
			{
			JsonObject ob = root.getAsJsonObject();
			final Interval interval = new Interval(
					ob.get("chrom").getAsString(),
					ob.get("start").getAsInt(),
					ob.get("end").getAsInt()
					);
			return Collections.singletonList(new IntervalSplit(interval));
			}
		else if(root.isJsonArray())
			{
			final JsonArray array=root.getAsJsonArray();
			if(array.size()==0) throw new IOException("Empty array");
			for(int i=0; i< array.size();++i) {
				JsonElement item = array.get(i);
				if(item.isJsonPrimitive())
					{
					splits.add(parseRefSplitFromStr(item.getAsString()));
					}
				else if(item.isJsonObject())
					{
					JsonObject ob = item.getAsJsonObject();
					final Interval interval = new Interval(
							ob.get("chrom").getAsString(),
							ob.get("start").getAsInt(),
							ob.get("end").getAsInt()
							);
					splits.add(new IntervalSplit(interval));
					}
				else
					{
					 throw new IOException("Illegal array item "+item);
					}
				}
			return splits;
			}
		else if(root.isJsonPrimitive())
			{
			final String str=root.getAsString();
			if(str.endsWith(".bed"))
				{
				splits.addAll(Files.readAllLines(Paths.get(str)).stream().
					filter(L->!L.startsWith("browser")).
					filter(L->!L.startsWith("track")).
					filter(L->!L.startsWith("#")).
					filter(L->!L.trim().isEmpty()).
					map(L->L.split("[\t]")).
					map(T->new Interval(T[0],1+Integer.parseInt(T[1]),Integer.parseInt(T[2]))).
					map(I->new IntervalSplit(I)).
					collect(Collectors.toList())
					);
				return splits;
				}
			else if(str.endsWith(".json"))
				{
				splits.addAll(parseRefSplitList(readJsonFile(new File(str))));
				return splits;
				}
			else if(str.endsWith(".intervals"))//gatk style
				{
				splits.addAll(Files.readAllLines(Paths.get(str)).stream().
						filter(L->!L.startsWith("#")).
						filter(L->!L.trim().isEmpty()).
						map(L-> parseRefSplitFromStr(L)).
						collect(Collectors.toList())
						);
				return splits;
				}
			else if(str.equals("genome_wide"))
				{
				return Collections.singletonList(new NoSplit());
				}
			else if(str.equals("by_chromosome"))
				{
				for(int i=1;i<=22;++i) splits.add(new ContigSplit(String.valueOf(i)));
				splits.add(new ContigSplit("X"));
				splits.add(new ContigSplit("Y"));
				return splits;
				}
			else 
				{
				return Collections.singletonList(parseRefSplitFromStr(str));
				}
			}
		else
			{
			throw new IOException("illegal refSplit "+root);
			}
		}
	
	private abstract class HasAttributes
		{
		private final Map<PropertyKey,String> _att=new HashMap<>();
		
		protected HasAttributes(final JsonElement root) throws IOException
			{
			if(root!=null && root.isJsonObject())
				{
				for(final PropertyKey prop:NgsWorkflow.name2propertyKey.values())
					{
					if(!root.getAsJsonObject().has(prop.getKey())) continue;
					this._att.put(prop, root.getAsJsonObject().get(prop.getKey()).getAsString());
					}
				}
			}
		
		public abstract HasAttributes getParent();
		public Map<PropertyKey,String> getAttributes() { return _att;}
		
		private boolean booleanValue(final PropertyKey key) {
			String s=getAttribute(key,"false");
			if(s!=null && s.equals("true")) return true;
			if(s!=null && s.equals("false")) return false;
			throw new RuntimeException("bad boolean value for "+key+" : \""+s+"\"");
			}
		public boolean isTrue(final PropertyKey key) {
			return booleanValue(key);
			}
		public boolean isFalse(final PropertyKey key) {
			return !booleanValue(key);
			}
		
		
		public String getAttribute(final PropertyKey key,final String def)
			{
			if(key==null) throw new NullPointerException("key is null");
			HasAttributes curr=this;
			while(curr!=null)
				{
				if(curr.getAttributes().containsKey(key))
					{
					return curr.getAttributes().get(key);
					}
				curr=curr.getParent();
				}
			if(def!=null) return def;
			if(key.getDefaultValue()!=null) return key.getDefaultValue();
			throw new RuntimeException("key "+key+" is not defined for "+this.toString());
			}
		
		public String getAttribute(final PropertyKey key)
			{
			return getAttribute(key,null);
			}
		
		public boolean isAttributeSet(final PropertyKey key,final Boolean def)
			{
			final String s= getAttribute(key,(def==null?null:def.booleanValue()?"true":"false"));
			if(s==null) throw new RuntimeException("key "+key+" is not defined for "+this.toString());
			if(s.equals("1") || s.equals("true") || s.equals("yes") || s.equals("T") || s.equals("Y")) return true;
			if(s.equals("0") || s.equals("false") || s.equals("no") || s.equals("F") || s.equals("N")) return false;
			throw new RuntimeException("bad boolean value for "+key+" : "+s);
			}
		
		public boolean isAttributeSet(final PropertyKey key)
			{
			return isAttributeSet(key,null);
			}
		
		public String getTmpPrefixToken()
			{
			return getAttribute(PROP_TMP_PREFIX_TOKEN);
			}
		public String getFilePrefix()
			{
			return getAttribute(PROP_FILES_PREFIX_TOKEN);
			}
		public String getTmpPrefix()
			{
			return getTmpPrefixToken()+getFilePrefix();
			}
		@Override
		public abstract String toString();
		
		protected String rulePrefix()
			{
			return "\tmkdir -p $(dir $@) && [[ ! -f \"${OUT}/STOP\" ]]  ";
			}

		}
	
	
	/** describe a NGS project */
	private class Project extends HasAttributes
		{
		private final String name;
		private final String description;
		private final List<Sample> _samples=new ArrayList<>();
		private final Optional<Capture> capture;
		private final Pedigree pedigree;
		private final Set<String> vcfToReannotate;
		
		Project(final JsonElement root) throws IOException
			{
			super(root);
			
			if(NgsWorkflow.this.reannotatevcfs)
				{
				this.vcfToReannotate=new LinkedHashSet<>();
				}
			else
				{
				this.vcfToReannotate=Collections.emptySet();
				}
			
			if(!root.isJsonObject()) throw new IOException("project json is not object");
			final JsonObject json = root.getAsJsonObject();
			
			if(!json.has("name")) {
				throw new IOException("project name missing");
			}
			this.name = json.get("name").getAsString();
			if(!json.has("description")) {
				this.description=this.name;
			} else
			{
				this.description = json.get("description").getAsString();
			}
			
			if(json.has("capture")) {
				final Capture c=new Capture(this,json.get("capture"));
				capture = Optional.of(c);
				}
			else
				{
				this.capture =Optional.empty();
				}
			
			if(json.has("pedigree")) {
				this.pedigree = new Pedigree(this,json.get("pedigree"));
				}
			else
				{
				this.pedigree = new Pedigree(this,new JsonObject());
				}
			
			final Set<String> seen=new HashSet<>(); 
			if(json.has("samples")) {
				for(final JsonElement samplejson :json.get("samples").getAsJsonArray())
					{
					Sample sample =new Sample(this,samplejson.getAsJsonObject());
					this._samples.add(sample);
					if(seen.contains(sample.getName()))
						{
						throw new IOException("Sample "+sample+" defined twice");
						}
					seen.add(sample.getName());
					}
				}
			else if(NgsWorkflow.this.reannotatevcfs && json.has("vcfs")) 
				{
				for(final JsonElement vcfjson :json.get("vcfs").getAsJsonArray())
					{
					if(vcfjson.isJsonPrimitive())
						{
						this.vcfToReannotate.add(vcfjson.getAsString());
						}
					else
						{
						throw new IOException("Not a string "+vcfjson);
						}
					}
				}
			else
				{
				LOG.warn("No 'samples' under project");
				}
			}
		
		public String getName() {
			return name;
			}
		
		public String getDescription() {
			return description;
		}
		
		public Pedigree getPedigree() {
			return pedigree;
		}
		
		public List<Sample> getSamples() { return this._samples;}
		@Override  public HasAttributes getParent() { return null;}
		
		public  String getOutputDirectory() {
			return getAttribute(PROP_PROJECT_OUTPUT_DIRECTORY);
			}
		public  String getSamplesDirectory() {
			return "${OUT}/Samples";
			}
		public  String getBedDirectory() {
			return "${OUT}/BED";
			}
		
		public  String getVcfDirectory() {
			return "${OUT}/VCF";
			}
		

		
		boolean hasCapture()
			{
			return this.capture.isPresent();
			}
		public Capture getCapture()
			{
			return this.capture.get();
			}
		
		/** when reannotation vcf, what will be the output filename */
		public String getReannoteVcfFilename( String s)
			{
			if(!NgsWorkflow.this.reannotatevcfs) throw new IllegalStateException();
			if(s==null || !(s.endsWith(".vcf") || s.endsWith(".vcf.gz"))) {
				throw new IllegalArgumentException("Bad extension vcf:"+s);
				}
			int slash= s.lastIndexOf('/');
			if(slash!=-1) s=s.substring(slash+1);
			if(s.matches("20[0-9][0-9][0-1][0-9][0-3][0-9].*"))
				{
				s=s.substring(8);
				}
			while(s.startsWith("_") || s.startsWith("."))
				{
				s=s.substring(1);
				}
			if(s.endsWith(".vcf"))
				{
				s=s.substring(0, s.length()-4);
				}
			else if(s.endsWith(".vcf.gz"))
				{
				s=s.substring(0, s.length()-7);
				}
			s+=".vcf.gz";
			return getOutputDirectory()+"/"+getFilePrefix()+s;
			}
		
		public  String getHapCallerAnnotationVcf() {
			return getVcfDirectory()+"/"+getFilePrefix()+"HCAnnotations.vcf.gz";
			}
		public  String getLumpyVcf() {
			return getVcfDirectory()+"/"+getFilePrefix()+"LumpyExpress.vcf.gz";
			}
		public  String getHapCallerGenotypedVcf() {
			return getVcfDirectory()+"/"+getTmpPrefix()+"HCGenotyped.vcf.gz";
			}
		public  String getSamtoolsRawVcf() {
			return getVcfDirectory()+"/"+getTmpPrefix()+"Samtools.vcf.gz";
			}
		
		public List<RefSplit> getHaplotypeCallerSplits() {
			final List<RefSplit> chroms=new ArrayList<>(25);
			for(int i=1;i<=22;++i) chroms.add(new ContigSplit(String.valueOf(i)));
			chroms.add(new ContigSplit("X"));
			chroms.add(new ContigSplit("Y"));
			return chroms;
			}
		
		public List<RefSplit> getSamtoolsCallerSplits() {
			return getHaplotypeCallerSplits();
			}
		
		

		String getFinalBamList() {
			return getSamplesDirectory()+"/"+getTmpPrefix()+"final.bam.list";
			}
		
		void bamList(final PrintWriter w) {
			w.println(getFinalBamList()+":"+this.getSamples().stream().map(S->S.getFinalBamBai()).collect(Collectors.joining(" \\\n\t")));
			w.print(rulePrefix()+ " && rm -f $@ $(addsuffix .tmp,$@) ");
			this.getSamples().stream().forEach(S->{
				w.print(" && echo \""+S.getFinalBam()+"\" >> $(addsuffix .tmp,$@) ");
				});
			w.println("&& mv --verbose $(addsuffix .tmp,$@) $@");
			}
		
		public List<AbstractCaller> getCallers()
			{
			final List<AbstractCaller> L = new ArrayList<>();
			L.add(new UnifiedGenotyperCaller());
			//L.add(new SamtoolsCaller());
			return L;
			}
		
		
		public void callVariants(PrintWriter out) {
			for(AbstractCaller C:getCallers()) C.print(out);
		}
		
		public void annotateVariants(PrintWriter out) {
			new VariantAnnotator(this.getHapCallerAnnotationVcf(),this.getHapCallerGenotypedVcf()).
				setPedigree( getPedigree().getPedFilename()).
				print(out);
		}
		
		void reannoteVcfs(PrintWriter out) {
			if(!NgsWorkflow.this.reannotatevcfs) throw new IllegalStateException();
			for(final String invcf: this.vcfToReannotate)
				{
				new VariantAnnotator(this.getReannoteVcfFilename(invcf),invcf).
					//setRemoveAnnotations().
					print(out);
				}
			}

		
		
		@Override
		public String toString() {
			return getName();
			}
		
		public String getNoResultContig() {
			return "22";
		}
		
		
		public class VariantAnnotator
			{
			private String target;
			private String dep;
			private String pedigree=null;
			private boolean removeAnnotations=false;
			public VariantAnnotator(String target,String dep)
				{
				this.target=target;
				this.dep=dep;
				}
			
			public VariantAnnotator setPedigree(final String ped)
				{
				if(!(ped==null || ped.trim().isEmpty()))
					{
					this.pedigree=ped;
					}
				return this;
				}
			public VariantAnnotator setRemoveAnnotations(boolean b) 
				{
				this.removeAnnotations=b;
				return this;
				}
			public VariantAnnotator setRemoveAnnotations() 
				{
				return setRemoveAnnotations(true);
				}
			
			void  print(final PrintWriter w)
				{
				w.print(this.target);
				w.print(":");
				w.print(this.dep);
				if(this.pedigree!=null)
					{
					w.print(" ");
					w.print(this.pedigree);
					}
				w.println();
				w.print(rulePrefix());
				w.print(" && rm -f  $(addsuffix .tmp1.vcf,$@)  $(addsuffix .tmp1.vcf.idx,$@) $(addsuffix .tmp2.vcf,$@)  $(addsuffix .tmp2.vcf.idx,$@) ");
				w.print(" && "+(this.dep.endsWith(".vcf")?"cat":"gunzip -c")+" $< | ");
				
				if( removeAnnotations ) {
					w.print(" $(call run_jvarkit,vcfstripannot) -x '"
							+ "INFO/DukeMapabilityUniqueness35bp,"
							+ "FILTER/DukeMapabilityUniqueness35LT1,"
							+ "FILTER/MapabilityConsensusExcludable,"
							+ "FILTER/MapabilityRegionsExcludable,"
							+ "FILTER/IN_EXAC,"
							+ "FILTER/IN_GNOMAD,"
							+ "INFO/ANN,"
							+ "INFO/CSQ,"
							+ "INFO/exac.AC,"
							+ "INFO/PossibleDeNovo"
							+ "' - | ");
					}
				
				w.print(" $(call run_jvarkit,vcfbigwig) "
						+ " -t ensembl2ucsc -T DukeMapabilityUniqueness35bp -B /commun/data/pubdb/ucsc/hg19/encodeDCC/wgEncodeDukeMapabilityUniqueness35bp.bigWig |");
				w.print(" $(call run_jvarkit,vcffilterjs) -F DukeMapabilityUniqueness35LT1 -e '!variant.hasAttribute(\"DukeMapabilityUniqueness35bp\") || variant.getAttributeAsDouble(\"DukeMapabilityUniqueness35bp\",100.0)>=1.0;' > $(addsuffix .tmp1.vcf,$@) ");
				//BED
				w.print(" && ${java.exe}   -Djava.io.tmpdir=$(dir $@)  -jar ${gatk.jar}  -T VariantFiltration -R $(REF) "
						+ " -o $(addsuffix .tmp2.vcf,$@) -L $(addsuffix .tmp1.vcf,$@) --variant $(addsuffix .tmp1.vcf,$@) "
						+" --maskName MapabilityConsensusExcludable --mask:BED /commun/data/pubdb/ucsc/hg19/encodeDCC/wgEncodeDacMapabilityConsensusExcludable_nochrprefix.bed.gz "
						+ " && mv --verbose   $(addsuffix .tmp2.vcf,$@)  $(addsuffix .tmp1.vcf,$@)"
						+ " && mv --verbose   $(addsuffix .tmp2.vcf.idx,$@)  $(addsuffix .tmp1.vcf.idx,$@)"
						);
				w.print(" && ${java.exe} -Djava.io.tmpdir=$(dir $@)  -jar ${gatk.jar}  -T VariantFiltration -R $(REF) "
						+ " -o $(addsuffix .tmp2.vcf,$@) -L $(addsuffix .tmp1.vcf,$@) --variant $(addsuffix .tmp1.vcf,$@) "
						+" --maskName MapabilityRegionsExcludable --mask:BED /commun/data/pubdb/ucsc/hg19/encodeDCC/wgEncodeDukeMapabilityRegionsExcludable_nochrprefix.bed.gz "
						+ " && mv --verbose   $(addsuffix .tmp2.vcf,$@)  $(addsuffix .tmp1.vcf,$@)"
						+ " && mv --verbose   $(addsuffix .tmp2.vcf.idx,$@)  $(addsuffix .tmp1.vcf.idx,$@)"
						);
				//EXAC
				w.print(" && ${java.exe}   -Djava.io.tmpdir=$(dir $@)  -jar ${gatk.jar}  -T VariantAnnotator -R $(REF) "
						+ " -o $(addsuffix .tmp2.vcf,$@) -L $(addsuffix .tmp1.vcf,$@) --variant $(addsuffix .tmp1.vcf,$@) "
						+" --resource:exac /commun/data/pubdb/broadinstitute.org/exac/1.0/ExAC.r1.sites.vcf.gz   --resourceAlleleConcordance "
						+" --expression exac.AC  ");
				if(this.pedigree!=null)
					{
					w.print(" --pedigree "+this.pedigree);
					w.print(" -A PossibleDeNovo ");
					}
				w.print(
						" -A HomopolymerRun -A TandemRepeatAnnotator  "
						+ " && mv --verbose   $(addsuffix .tmp2.vcf,$@)  $(addsuffix .tmp1.vcf,$@)"
						+ " && mv --verbose   $(addsuffix .tmp2.vcf.idx,$@)  $(addsuffix .tmp1.vcf.idx,$@)"
						);
				// filter exac
				w.print(" && $(call run_jvarkit,vcffilterjs) -F IN_EXAC -e '!variant.hasAttribute(\"exac.AC\")' $(addsuffix .tmp1.vcf,$@) |");
				//gnomad
				w.print(" $(call run_jvarkit,vcfgnomad) -ac -gf IN_GNOMAD -m '${gnomad.vcf.manifest}' |");
				//snpeff
				w.print(" ${java.exe}  -Djava.io.tmpdir=$(dir $@) -jar ${snpeff.jar} ann -c ${snpeff.config} GRCh37.75 -nodownload -noStats |  ");
				//vep
				w.print(" $(call vep78) |");
				w.print(" ${bgzip.exe} >  $(addsuffix .tmp2.vcf.gz,$@)  ");
				w.print(" && ${tabix.exe} -p vcf -f $(addsuffix .tmp2.vcf.gz,$@) "
					+ " && mv --verbose   $(addsuffix .tmp2.vcf.gz,$@) $@ "
					+ " && mv --verbose   $(addsuffix .tmp2.vcf.gz.tbi,$@)  $(addsuffix .tbi,$@)"
					);
				w.print(" && rm -f  $(addsuffix .tmp1.vcf,$@)  $(addsuffix .tmp1.vcf.idx,$@)  $(addsuffix .tmp2.vcf,$@)  $(addsuffix .tmp2.vcf.idx,$@) ");
				w.println();
				}
			}
		
		
		private class LumpyExpress
			{
			private List<LumpyInput> inputs=new ArrayList<>();
			Project getProject() { return Project.this;}
			
			LumpyExpress addBam(final String inputBam)
				{
				this.inputs.add(new LumpyInput(inputBam));
				return this;
				}
			
			private  class LumpyInput
				{
				final String inputBam;
				LumpyInput(final String inputBam) {
					this.inputBam = inputBam;
					}
				
				private String getDiscordantBam() 
					{
					return "$(addsuffix .tmp.discordant.bam,"+inputBam+")";
					}
				private String getSplitterBam() 
					{
					return "$(addsuffix .tmp.splitters.bam,"+inputBam+")";
					}
				
				void  print(final PrintWriter w)  {
					// Extract the discordant paired-end alignments.
					w.print(this.getDiscordantBam());
					w.print(":");
					w.print(this.inputBam);
					w.println();
					w.print(rulePrefix());
					w.print(" && ${samtools.exe} view -F 1294 -b -o $(addsuffix .tmp.bam,$@) $< ");
					w.print(" && mv --verbose $(addsuffix .tmp.bam,$@) $@ ");
					w.println();
					
					// Extract the split-read alignments
					w.print(this.getSplitterBam());
					w.print(":");
					w.print(this.inputBam);
					w.println();
					w.print(rulePrefix());
					w.print(" && ${samtools.exe} view -h  $< | ");
					w.print(" ${lumpy.dir}/scripts/extractSplitReads_BwaMem -i stdin |");
					w.print(" ${samtools.exe} view -Sb -o $(addsuffix .tmp.bam,$@) - "); 
					w.print(" && mv --verbose $(addsuffix .tmp.bam,$@) $@ ");
					w.println();
					}
				}
			void  print(final PrintWriter w)
				{
				if(this.inputs.isEmpty()) return;
				
				for(LumpyInput input:this.inputs) {
					input.print(w);
					}
				
				w.print(this.getProject().getLumpyVcf());
				w.print(":");
				w.println(this.inputs.stream().map(T->T.getDiscordantBam()+" "+T.getSplitterBam()+ " "+T.inputBam).collect(Collectors.joining(" ")));
				w.println();
				w.print(rulePrefix());
				w.print("&&  ${lumpyexpress.exe} -B ");
				w.print(this.inputs.stream().map(T->T.inputBam.trim()).collect(Collectors.joining(",")));
				w.print(" -S ");
				w.print(this.inputs.stream().map(T->T.getSplitterBam().trim()).collect(Collectors.joining(",")));
				w.print(" -D ");
				w.print(this.inputs.stream().map(T->T.getDiscordantBam().trim()).collect(Collectors.joining(",")));
				w.print(" -o $(addsuffix .tmp.vcf,$@)  ");
				w.print(" && $(java.exe)  -Djava.io.tmpdir=$(dir $@) -jar $(picard.jar) UpdateVcfSequenceDictionary I=$(addsuffix .tmp.vcf,$@)  O=$(addsuffix .tmp2.vcf,$@)   SEQUENCE_DICTIONARY=$(addsuffix .dict,$(basename ${REF})) && mv --verbose $(addsuffix .tmp2.vcf,$@) $(addsuffix .tmp.vcf,$@) ");
				w.print(" && $(java.exe)  -Djava.io.tmpdir=$(dir $@) -jar $(picard.jar) SortVcf I=$(addsuffix .tmp.vcf,$@)  O=$(addsuffix .tmp2.vcf,$@)  && mv --verbose $(addsuffix .tmp2.vcf,$@) $(addsuffix .tmp.vcf,$@) ");				
				w.print(" && ${bgzip.exe} -f $(addsuffix .tmp.vcf,$@)   ");
				w.print(" && ${tabix.exe} -p vcf -f $(addsuffix .tmp.vcf.gz,$@) "
					+ " && mv --verbose   $(addsuffix .tmp.vcf.gz,$@) $@ "
					+ " && mv --verbose   $(addsuffix .tmp.vcf.gz.tbi,$@)  $(addsuffix .tbi,$@)"
					);
				w.print(" && rm --verbose "+ this.inputs.stream().map(T->T.getDiscordantBam()+" "+T.getSplitterBam()).collect(Collectors.joining(" ")));
				w.print(" && touch  "+ this.inputs.stream().map(T->T.getDiscordantBam()+" "+T.getSplitterBam()).collect(Collectors.joining(" ")));
				w.print(" && sleep 2 && touch -c $@  $(addsuffix .tbi,$@)");
				w.println();
				}
			
			}

		void lumpyExpress(final PrintWriter out) {
			if(this.isFalse(PROP_USE_LUMPY_EXPRESS)) return;
			final LumpyExpress lumpy=new LumpyExpress();
			for(Sample sample:this.getSamples()) {
				lumpy.addBam(sample.getFinalBam());
			}
			lumpy.print(out);
		}
		
		
		private abstract class AbstractCaller
			{
			Project getProject() { return Project.this;}
			abstract List<? extends RefSplit> getCallSplits();
			abstract String getTargetVcfFilename();
			
			abstract void call(final PrintWriter w,RefSplit split);
			
			void  print(final PrintWriter w)  {
				final List<String> vcfParts=new ArrayList<>();
				final List<? extends RefSplit> refSplits = this.getCallSplits();
				if(refSplits.isEmpty()) throw new IllegalArgumentException();
				
				for(final RefSplit split: refSplits)
					{
					final String vcfPart;
					switch(split.getType()) {
						case WHOLE_GENOME:  vcfPart= getTargetVcfFilename();break;
						default:  vcfPart= "$(addsuffix "+split.getToken()+".vcf.gz,"+getTargetVcfFilename()+")";
						}
					
					vcfParts.add(vcfPart);
					
					w.append(vcfPart).append(":").append(getFinalBamList());
					w.append(" ").append(getPedigree().getPedFilename());
					if( getProject().hasCapture())
		            	{
						w.append(" ").append(getCapture().getExtendedFilename());
		            	}
					w.append("\n");
					w.append(rulePrefix()+" && ");
					
					if( getProject().hasCapture())
		            	{
						switch(split.getType())
							{
							case WHOLE_GENOME: break;
							case INTERVAL:{
									final IntervalSplit tmp= IntervalSplit.class.cast(split);
									w.append(" ${bedtools.exe} intersect -a ").append(getCapture().getExtendedFilename()).
										append(" -b <(echo -e '").append(tmp.getInterval().getContig()).append("\\t").
										append(String.valueOf(tmp.getInterval().getStart())).append("\\t").
										append(String.valueOf(tmp.getInterval().getEnd())).append("') | ").
										append(" awk -F '\t' 'BEGIN{N=0;}{print;N++;}END{if(N==0) printf(\""+getNoResultContig()+"\\t0\\t1\\n\");}' > $(addsuffix .bed,$@) && "); 
									break;
									}
							case WHOLE_CONTIG:
									{
									final ContigSplit tmp= ContigSplit.class.cast(split);
									w.append(" awk -F '\t' 'BEGIN{N=0;}{if($$1==\""+tmp.getContig()+"\") {print;N++;}}END{if(N==0) printf(\""+getNoResultContig()+"\\t0\\t1\\n\");}' "+getCapture().getExtendedFilename()+" > $(addsuffix .bed,$@) && ");
									break;
									}
							default: throw new IllegalStateException();
							}
		            	}
					this.call(w,split);
					w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz,$@)\" \"$@\" ");
					w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz.tbi,$@)\" \"$(addsuffix .tbi,$@)\" ");
					
					if( getProject().hasCapture())
		            	{
						w.append(" && rm --verbose \"$(addsuffix .bed,$@)\" ");
		            	}
					w.append("\n");
					}
				
				if(vcfParts.size()>1) {
					w.append(getTargetVcfFilename()).append(":");
					for(final String vcfPart:vcfParts) {
						w.append(" \\\n\t").append(vcfPart);
						}
					w.append("\n");
					w.append(rulePrefix()+" && rm -f $(addsuffix .list,$@) ");
					int vcfn=0;
					for(final String vcfPart:vcfParts) {
						if(vcfn%30==0)
							{
							w.append("\n\t");
							}
						else
							{
							w.append(" && ");
							}
						w.append("echo '");
						w.append(vcfPart);
						w.append("' >>  $(addsuffix .list,$@) ");
						vcfn++;
						}
					
					w.append("\n\t${java.exe}   -Djava.io.tmpdir=$(dir $@)  -jar ${gatk.jar}  -T CombineVariants -R $(REF) "
							+ " -o $(addsuffix .tmp.vcf.gz,$@) -genotypeMergeOptions UNSORTED "
							+" --variant  $(addsuffix .list,$@) "
							);
					
					w.append(" && rm --verbose $(addsuffix .list,$@)");
					w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz,$@)\" \"$@\" ");
					w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz.tbi,$@)\" \"$(addsuffix .tbi,$@)\" ");
					w.append("\n");
					}
				}
			}
	
	
	private class UnifiedGenotyperCaller extends AbstractCaller
		{	
		@Override
		List<? extends RefSplit> getCallSplits() {
			return getIntervalSplitListGapForGrch37();
			}
		@Override String getTargetVcfFilename() { return  getProject().getHapCallerGenotypedVcf();}
		
		@Override void call(final PrintWriter w,final RefSplit split)
			{
			w.append(" $(java.exe)  -XX:ParallelGCThreads=5 -Xmx2g  -Djava.io.tmpdir=$(dir $@) -jar $(gatk.jar) -T HaplotypeCaller ");
			w.append("	-R $(REF) ");
			w.append("	--validation_strictness LENIENT ");
			w.append("	-I $< -o \"$(addsuffix .tmp.vcf.gz,$@)\" ");
			//w.print("	--num_cpu_threads_per_data_thread "+  this.getIntProperty("base-recalibrator-nct",1));
			w.append("	-l INFO ");
			w.append("	-nct ").append(getAttribute(PROP_GATK_HAPCALLER_NCT));
			
			w.append("	--dbsnp \"$(gatk.bundle.dbsnp.vcf)\" ");
			w.append(" $(foreach A, StrandAlleleCountsBySample  PossibleDeNovo AS_FisherStrand AlleleBalance AlleleBalanceBySample "
					/* BaseCountsBySample fait planter */
					+ " GCContent ClippingRankSumTest , --annotation ${A} ) ");
			w.append(" --pedigree ").append(getProject().getPedigree().getPedFilename());
			if( getProject().hasCapture())
	            {
				switch(split.getType())
					{
					case WHOLE_GENOME: 	w.append(" -L:"+getProject().getCapture().getName()+",BED "+getCapture().getExtendedFilename());
					case INTERVAL: //through...
					case WHOLE_CONTIG: w.append(" -L:BED \"$(addsuffix .bed,$@)\" ");break;
					default: throw new IllegalStateException();
					}
	            }
			else
				{
				switch(split.getType())
					{
					case WHOLE_GENOME: /* nothing */ break;
					case INTERVAL: {
						IntervalSplit tmp= IntervalSplit.class.cast(split);
						w.append(" -L \""+tmp.getInterval().getContig()+":"+ tmp.getInterval().getStart()+"-"+tmp.getInterval().getEnd()+"\" ");
						break;
						}
					case WHOLE_CONTIG: w.append(" -L ").append(ContigSplit.class.cast(split).getContig());break;
					default: throw new IllegalStateException();
					}					
				}
			}

		}

	private class SamtoolsCaller extends AbstractCaller
		{	
		@Override
		List<RefSplit> getCallSplits() {
			return getProject().getSamtoolsCallerSplits();
			}

		@Override String getTargetVcfFilename() { return  getProject().getSamtoolsRawVcf();}
		
		@Override void call(final PrintWriter w,final RefSplit split)
			{
			w.append(" ${samtools.exe} mpileup --uncompressed --BCF --fasta-ref $(REF) --bam-list $< ");
			
			if( getProject().hasCapture())
	            {
				switch(split.getType())
					{
					case WHOLE_GENOME: 	w.append(" --positions "+getCapture().getExtendedFilename());
					case INTERVAL: //through...
					case WHOLE_CONTIG: w.append(" --positions $(addsuffix .bed,$@) ");break;
					default: throw new IllegalStateException();
					}
	            }
			else
				{
				switch(split.getType())
					{
					case WHOLE_GENOME: /* nothing */ break;
					case INTERVAL: {
						final IntervalSplit tmp= IntervalSplit.class.cast(split);
						w.append(" --region \""+tmp.getInterval().getContig()+":"+ tmp.getInterval().getStart()+"-"+tmp.getInterval().getEnd()+"\" ");
						break;
						}
					case WHOLE_CONTIG: w.append(" --region ").append(ContigSplit.class.cast(split).getContig());break;
					default: throw new IllegalStateException();
					}					
				}
			
			
			w.append(" | ${bcftools.exe} call -vmO z -o \"$(addsuffix .tmp.vcf.gz,$@)\" ");
			w.append(" --samples-file ").append(getProject().getPedigree().getPedFilename());
			w.append(" && ${tabix.exe} -f -p vcf \"$(addsuffix .tmp.vcf.gz,$@)\" ");
			}
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		}
	
	/** describe a Pedigree capture */
	private class Pedigree extends HasAttributes
		{
		private final Project project;
		private final String pedFilename;
		private final boolean providedByUser;
		Pedigree( final Project project,final JsonElement root) throws IOException
			{
			super(root);
			this.project=project;
			if(root!=null && root.isJsonObject() && root.getAsJsonObject().has("ped"))
				{
				final JsonObject json=root.getAsJsonObject();
				this.pedFilename=json.get("ped").getAsString();
				this.providedByUser = true;
				}
			else if(root!=null && root.isJsonPrimitive())
				{
				this.pedFilename=root.getAsString();
				this.providedByUser = true;
				}
			else
				{
				this.providedByUser = false;
				this.pedFilename = project.getOutputDirectory()+"/PED/"+getTmpPrefix()+"pedigree.ped";
				}
			}
		public boolean isProvidedByUser() {
			return providedByUser;
			}
		public String getPedFilename() {
			return pedFilename;
			}
		public Project getProject() {
			return project;
			}
		@Override
		public Project getParent() {
			return getProject();
			}
		void build(final PrintWriter w)
			{
			if(isProvidedByUser()) return;
			w.println(getPedFilename()+":");
			w.println(rulePrefix()+" && rm --verbose -f $@ $(addsuffix .tmp.ped,$@)");
			for(final Sample sample: getProject().getSamples())
				{
				w.print("\techo '");
				w.print(String.join("\t",sample.getFamilyId(),sample.getName(),sample.getFatherId(),sample.getMotherId(),sample.getSex(),sample.getStatus()));
				w.println("' >> $(addsuffix .tmp.ped,$@)");
				}
			w.println("\tmv  --verbose $(addsuffix .tmp.ped,$@) $@");
			}
		
		@Override
		public String toString() {
			return getPedFilename();
			}
		}
	
	/** describe a BED capture */
	private class Capture extends HasAttributes
		{
		private final Project project;
		private Integer extend=null;
		private String name="capture";
		private final String bed;
		
		Capture( final Project project,final JsonElement root) throws IOException
			{
			super(root);
			this.project=project;
			if(root.isJsonObject())
				{
				final JsonObject json=root.getAsJsonObject();
				this.bed=json.get("bed").getAsString();
				if(json.has("extend"))
					{
					this.extend=json.get("extend").getAsInt();
					}
				if(json.has("name"))
					{
					this.name=json.get("name").getAsString();
					}
				}
			else if(root.isJsonPrimitive())
				{
				bed=root.getAsString();
				}
			else
				{
				throw new IOException("bad capture");
				}
			}
		
		public String getBedFilename()
			{
			return this.bed;
			}
		public String getName() {
			return name;
			}
		public String getNormalizedFilename()
			{
			return getProject().getBedDirectory()+"/"+getName()+".normalized.bed";
			}
		
		public String getExtendedFilename()
			{
			return getProject().getBedDirectory()+"/"+getName()+".extended"+getExtend()+".bed";
			}
		
		public int getExtend()
			{
			if(this.extend!=null) return this.extend;
			return Integer.parseInt(getAttribute(PROP_CAPTURE_EXTEND));
			}
		
		public Project getProject() {
			return project;
			}
		@Override
		public HasAttributes getParent() {
			return getProject();
			}
		
		
		void prepareCapture(final PrintWriter w)
		{
			w.println(getExtendedFilename()+": "+getNormalizedFilename()+" $(addsuffix .fai,$(REF))");
			w.print("\tmkdir -p $(dir $@) &&");
			w.print(" $(bedtools.exe) slop -b "+getExtend()+" -g $(addsuffix .fai,$(REF)) -i $< |");
			w.print(" LC_ALL=C sort -t '\t' -k1,1 -k2,2n -k3,3n |");
			w.print(" $(bedtools.exe) merge  -d  "+getExtend()+" |");
			w.print(" $(bedtools.exe) sort -faidx $(addsuffix .fai,$(REF)) > $(addsuffix .tmp.bed,$@) && ");
			w.print(" mv --verbose $(addsuffix .tmp.bed,$@) $@");
			w.println();

			w.println(getNormalizedFilename()+": "+getBedFilename()+" $(addsuffix .fai,$(REF))");
			w.print("\tmkdir -p $(dir $@) && ");
			w.print("	tr -d '\\r' < $< | grep -v '^browser'  | grep -v '^track' | cut -f1,2,3 | sed 's/^chr//' | LC_ALL=C sort -t '\t' -k1,1 -k2,2n -k3,3n | ");
			w.print("	$(bedtools.exe) merge |");
			w.print("   $(bedtools.exe) sort -faidx $(addsuffix .fai,$(REF)) > $(addsuffix .tmp.bed,$@) && ");
			w.print("	mv --verbose $(addsuffix .tmp.bed,$@) $@");
			w.println();

			}
		@Override
		public String toString() {
			return getBedFilename();
			}
		}
	
	private class Sample extends HasAttributes
		{
		private final Project project;
		private final String name;
		private final List<PairedFastq> pairs;
		private final String finalBam;
		private final String urer_g_vcf;
		
		private boolean isEmpty(File f) throws IOException
			{
			InputStream in=null;
			try {
				in = new FileInputStream(f);
				if(f.getName().endsWith(".gz")) in=new GZIPInputStream(in);
				int c=in.read();
				in.close();
				in=null;
				return c==-1;
			} catch (Exception e) {
				throw new RuntimeIOException(e);
			} finally {
				CloserUtil.close(in);
			}
			}
		
		Sample( final Project project,final JsonObject root) throws IOException
			{
			super(root);
			if(!root.isJsonObject()) throw new IOException("sample json is not object");
			final JsonObject json = root.getAsJsonObject();
	
			this.project = project;
			
			if(!json.has("name")) throw new IOException("@name missing in sample");
			this.name= json.get("name").getAsString();
			if(json.has("bam")) {
				this.pairs = Collections.emptyList();
				this.finalBam=json.get("bam").getAsString();
				}
			else if(json.has("fastqs"))
				{
				this.finalBam=null;
				this.pairs = new ArrayList<>();
				final JsonElement fastqs = json.get("fastqs");
				if(fastqs.isJsonPrimitive()) {
					final File dir=new File(fastqs.getAsString());
					if(!(dir.exists() && dir.isDirectory())) throw new IOException("not a directory "+dir );
					final File inside[]=dir.listFiles(F->F.exists() && F.getName().contains("_R1_") && F.getName().endsWith(".fastq.gz"));
					if(inside==null || inside.length==0) throw new IOException("No file under "+dir);
					for(final File insideFileR1:inside)
						{
						if(isEmpty(insideFileR1)) {
							LOG.warn(insideFileR1.getPath()+" is empty!");
							continue;
						}
						final JsonArray twofastqs=new JsonArray();
						twofastqs.add(new JsonPrimitive(insideFileR1.getAbsolutePath()));
						final File insideFileR2=new File(dir,insideFileR1.getName().replaceAll("_R1_", "_R2_"));
						if(!insideFileR2.exists()) throw new IOException("Cannot find "+insideFileR2);
						twofastqs.add(new JsonPrimitive(insideFileR2.getAbsolutePath()));
						final PairedFastq pair=new PairedFastq(this, this.pairs.size(), twofastqs);
						this.pairs.add(pair);
						}
					}
				else if(fastqs.isJsonArray())
					{
					for(final JsonElement pairjson : fastqs.getAsJsonArray())
						{
						final PairedFastq pair =new PairedFastq(this,this.pairs.size(),pairjson);
						this.pairs.add(pair);
						}
					}
				else
					{
					throw new IOException("bad fastq");
					}
				}
			else
				{
				this.finalBam=null;
				this.pairs = Collections.emptyList();
				}
			
			if(json.has("gvcf"))
				{
				this.urer_g_vcf=json.get("gvcf").getAsString();
				}
			else
				{
				this.urer_g_vcf=null;
				}
			}
		
		public String getFamilyId() { return "PEDIGREE";}
		public String getFatherId() { return "0";}
		public String getMotherId() { return "0";}
		public String getSex() { return "0";}
		public String getStatus() { return "9";}
		
		
		public Project getProject() { return this.project;}
		public List<PairedFastq> getPairs() { return this.pairs;}
		public String getName() { return this.name;}
		public boolean isHaloplex() { return isAttributeSet(PROP_PROJECT_IS_HALOPLEX, false);}
		public boolean isBamAlreadyProvided() { return this.finalBam!=null;}
		
		@Override public Project getParent() { return getProject();}
		
		public String getDirectory() {
			return getProject().getSamplesDirectory()+"/"+getName();
			}
		
		public String getMergedBam() {
			return getDirectory()+"/"+getTmpPrefix()+getName()+".merged.bam";
			}
		public String getMergedBamBai() {
			return "$(call  picardbai,"+getMergedBam()+")";
			}
		
		public boolean isMarkupDisabled() {
			return this.finalBam!=null || isHaloplex() || isAttributeSet(PROP_DISABLE_MARKDUP, false);
			}
		
		public String getMarkdupBam() {
			if(isMarkupDisabled()) {
				return getMergedBam();
				}
			else
				{
				return getDirectory()+"/"+getTmpPrefix()+getName()+".markdup.bam";
				}
			}
		public String getMarkdupBamBai() {
			return "$(call picardbai,"+getMarkdupBam()+")";
			}

		public boolean isRealignDisabled() {
			return  isBamAlreadyProvided() || isAttributeSet(PROP_DISABLE_REALIGN, false);
			}

		
		public String getRealignBam() {
			if(isRealignDisabled()) {
				return getMarkdupBam();
				}
			else
				{
				return getDirectory()+"/"+getTmpPrefix()+getName()+".realign.bam";
				}
			}
		public String getRealignBamBai() {
			return "$(call picardbai,"+getRealignBam()+")";
			}

		
		public boolean isRecalibrationDisabled()
			{
			return isBamAlreadyProvided() || isAttributeSet(PROP_DISABLE_RECALIBRATION, false);
			}
		
		
		public String getRecalBam() {
			if(isRecalibrationDisabled()) {
				return getRealignBam();
				}
			else
				{
				return getDirectory()+"/"+getTmpPrefix()+getName()+".recal.bam";
				}
			}
		public String getRecalBamBai() {
			return "$(call picardbai,"+getRecalBam()+")";
			}
		
		public String getFinalBam() {
			if(isBamAlreadyProvided()) {
				return this.finalBam;
				}
			else
				{
				return getDirectory()+"/"+getFilePrefix()+getName()+".final.bam";
				}
			}
		public String getFinalBamBai() {
			return "$(call picardbai,"+getFinalBam()+")";
			}
		
		public String mergeSortedBams()
			{
			if(isBamAlreadyProvided()) return "";
			return new StringBuilder().
				append(getMergedBamBai()+" : "+getMergedBam()).append("\n\ttouch -c $@\n").
				
				append(getMergedBam()+" : "+getPairs().stream().map(P->P.getSortedFilename()).collect(Collectors.joining(" "))+"\n").
				append(rulePrefix()+" && ").
				append("$(java.exe)  -Djava.io.tmpdir=$(dir $@) -jar \"$(picard.jar)\" MergeSamFiles O=$(addsuffix .tmp.bam,$@) ").
				append(" SO=coordinate AS=true CREATE_INDEX=true  COMPRESSION_LEVEL=").append(getAttribute(PROP_DEFAULT_COMPRESSION_LEVEL)).
				append(" VALIDATION_STRINGENCY=SILENT USE_THREADING=true VERBOSITY=INFO TMP_DIR=$(dir $@) COMMENT=\"Merged from $(words $^) files.\" ").
				append("$(foreach B,$(filter %.bam,$^), I=$(B) ) && ").
				append("mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\" ").
				append("\n").
				toString();
			}
		
		public String markDuplicates()
			{
			if(this.isMarkupDisabled()) return "";
			
			return new StringBuilder().
					append(getMarkdupBamBai()+" : "+getMarkdupBam()+"\n").
					append("\ttouch -c $@").
					append("\n").
					append(getMarkdupBam()+" : "+getMergedBam()+"\n").
					append(rulePrefix()+" && ").
					append("$(java.exe) -Xmx3g -Djava.io.tmpdir=$(dir $@) -jar \"$(picard.jar)\" MarkDuplicates I=$< O=$(addsuffix .tmp.bam,$@) M=$(addsuffix .metrics,$@) ").
					append(" REMOVE_DUPLICATES=").
					append(isAttributeSet(PROP_PICARD_MARKDUP_REMOVES_DUPLICATES)?"true":"false").
					append(" AS=true CREATE_INDEX=true  COMPRESSION_LEVEL=").
					append(getAttribute(PROP_DEFAULT_COMPRESSION_LEVEL)).
					append(" VALIDATION_STRINGENCY=SILENT VERBOSITY=INFO TMP_DIR=$(dir $@)  ").
					append(" && mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\" ").
					append(" && mv --verbose \"$(addsuffix .tmp.bai,$@)\" \""+getMarkdupBamBai()+"\"").
					append("\n").
					toString();
			}
		
		
		public String finalBam()
			{
			if(this.isBamAlreadyProvided()) return "";
			StringBuilder sb=new StringBuilder().
				append(getFinalBamBai()+":"+ this.getFinalBam()).
				append("\n").
				append(rulePrefix()).
				append(" && ${samtools.exe} index $< $@").
				append("\n")
				;
			sb.append(getFinalBam()+":"+ this.getRecalBam()).
				append("\n").
				append(rulePrefix()).
				append(" && cp --verbose $< $@").
				append("\n")
				;
			
			return sb.toString();
			}
		
		public String recalibrateBam()
		{
			
			if(this.isRecalibrationDisabled()) return "";

			final StringWriter sw=new StringWriter();
			final PrintWriter w= new PrintWriter(sw);

			w.println( getRecalBamBai() + " : " + this.getRecalBam());
			w.println("\ttouch -c $@");
			w.println();
		
		

			w.print(this.getRecalBam()+"  : "+this.getRealignBam()+" "+this.getRealignBamBai());
			if( this.getProject().hasCapture())
				{
				w.print(" " +getProject().getCapture().getExtendedFilename());
				}
			w.println();
			/* first call BaseRecalibrator */
			w.print("\tmkdir -p $(dir $@) && ");
			w.print(" $(java.exe)  -XX:ParallelGCThreads=5 -Xmx2g  -Djava.io.tmpdir=$(dir $@) -jar $(gatk.jar) -T BaseRecalibrator ");
			w.print("	-R $(REF) ");
			w.print("	--validation_strictness LENIENT ");
			w.print("	-I $< ");
			//w.print("	--num_cpu_threads_per_data_thread "+  this.getIntProperty("base-recalibrator-nct",1));
			w.print("	-l INFO ");
			w.print("	-o  $(addsuffix .grp,$@) ");
			w.print("	-knownSites:dbsnp138,VCF \"$(gatk.bundle.dbsnp.vcf)\" ");
			if( this.getProject().hasCapture())
                {
				w.print(" -L:"+getProject().getCapture().getName()+",BED "+getProject().getCapture().getExtendedFilename());
                }
			w.print("	-cov QualityScoreCovariate ");
			w.print("	-cov CycleCovariate ");
			w.print("	-cov ContextCovariate ");
			
			/* then call PrintReads */
			w.print(" && $(java.exe)  -XX:ParallelGCThreads=5 -Xmx4g  -Djava.io.tmpdir=$(dir $@) -jar $(gatk.jar) -T PrintReads ");
			w.print("	-R $(REF) ");
			w.print(" --bam_compression "+getAttribute(PROP_DEFAULT_COMPRESSION_LEVEL));
			w.print("	-BQSR  $(addsuffix .grp,$@) ");
			
			
			w.print(" --disable_indel_quals ");
				
			w.print("	-I $< ");
			w.print("	-o $(addsuffix .tmp.bam,$@) ");
			w.print("	--validation_strictness LENIENT ");
			w.print("	-l INFO && ");
			w.print("    mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\" && ");
			w.print("    mv --verbose \"$(call picardbai,$(addsuffix .tmp.bam,$@))\" \""+getRecalBamBai()+"\" && ");
			w.print("    sleep 2 && touch -c \""+getRecalBamBai()+"\" && ");
			w.print("  rm --verbose -f  $(addsuffix .grp,$@)");
			w.println();
			w.println();
			
			return sw.toString();
		}
		
		public String realignAroundIndels()
			{
			if(isRealignDisabled()) return "";
			
			StringWriter sw=new StringWriter();
			PrintWriter w= new PrintWriter(sw);
			
			w.append(getRealignBamBai()+" : "+getRealignBam()+"\n").
				append("\ttouch -c $@").
				append("\n")
				;
			
			/* first call RealignerTargetCreator */
			w.print(this.getRealignBam());
			w.print(" : ");
			w.print(this.getMarkdupBam());
			if( this.getProject().hasCapture())
				{
				w.print(" " +getProject().getCapture().getExtendedFilename());
				}
			w.println();
			w.print("\tmkdir -p $(dir $@ ) && ");
			w.print("$(java.exe)  -XX:ParallelGCThreads=2  -Xmx2g  -Djava.io.tmpdir=$(dir $@) -jar $(gatk.jar) -T RealignerTargetCreator ");
			w.print(" -R $(REF)  ");
			w.print(" --num_threads 1");
			if( this.getProject().hasCapture())
				{
				w.print(" -L:"+getProject().getCapture().getName()+",BED "+getProject().getCapture().getExtendedFilename());
				}
			w.print(" -I $(filter %.bam,$^) ");
			w.print(" -o $(addsuffix .intervals,$@) ");
			w.print(" --known:onekindels,VCF \"$(gatk.bundle.onek.indels.vcf)\" ");
			w.print(" --known:millsindels,VCF \"$(gatk.bundle.mills.indels.vcf)\"  ");
			w.print(" -S SILENT ");
			/* then call IndelRealigner */
			w.print("  && $(java.exe)  -XX:ParallelGCThreads=5 -Xmx2g  -Djava.io.tmpdir=$(dir $@) -jar $(gatk.jar) -T IndelRealigner ");
			w.print("  -R $(REF) ");
			w.print(" --bam_compression "+ getAttribute(PROP_DEFAULT_COMPRESSION_LEVEL));
			w.print("  -I $(filter %.bam,$^) "); 
			w.print("  -o $(addsuffix .tmp.bam,$@)  ");
			w.print("  -targetIntervals  \"$(addsuffix .intervals,$@)\" ");
			w.print("  -known:onekindels,VCF \"$(gatk.bundle.onek.indels.vcf)\"   ");
			w.print("  -known:millsindels,VCF \"$(gatk.bundle.mills.indels.vcf)\"  ");
			w.print("  -S SILENT &&  ");
			w.print("mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\" &&  ");
			w.print("mv --verbose  \"$(call picardbai,$(addsuffix .tmp.bam,$@))\" \""+ getRealignBamBai() +"\" &&  ");
			w.print("rm --verbose -f \"$(addsuffix .intervals,$@)\" ");
			w.println();
			w.println();

			return sw.toString();
			}
		
		@Override
		public String toString() {
			return getName();
			}
		}
	
	
	
	private class PairedFastq  extends HasAttributes
		{
		private final List<Fastq> fastqs;
		private final int indexInSample;
		private final Sample _sample;
	
		PairedFastq( final Sample sample,int indexInSample,final JsonElement root) throws IOException
			{
			super(root);
			this.indexInSample=indexInSample;
			this.fastqs = new ArrayList<>();
			this._sample=sample;
			if(root.isJsonObject())
				{
				final JsonObject json=root.getAsJsonObject();
				if(json.has("fastqs"))
					{
					for(final JsonElement fqjson :json.get("fastqs").getAsJsonArray())
						{		
						Fastq fq =new Fastq(this,fqjson);
						this.fastqs.add(fq);
						}
					}
				}
			else if(root.isJsonArray())
				{
				for(final JsonElement fqjson :root.getAsJsonArray())
					{		
					Fastq fq =new Fastq(this,fqjson);
					this.fastqs.add(fq);
					}
				}
			else if(root.isJsonPrimitive())
				{
				Fastq fq =new Fastq(this,root);
				this.fastqs.add(fq);
				}
			else throw new IOException("cannot decore paired for "+root);
			}
		
		public Sample getSample() { return this._sample;}
		public Project getProject() { return getSample().getProject();}
		
		@Override
		public final Sample getParent() { return getSample(); }
	
		public Fastq get(int i) { return this.fastqs.get(i);}
		public int getIndex() { return indexInSample;}
		
		
		public  String getDirectory() {
			return getSample().getDirectory()+"/p"+getIndex();
			}
		
		public String getSortedFilename()
		 	{
			if(getSample().getPairs().size()==1)
				{
				return getSample().getMergedBam();
				}
			else
				{
				return getDirectory()+"/"+getTmpPrefix()+getSample().getName()+".sorted.bam";
				}
			}
		private boolean isUsingCutadapt()
			{
			return isAttributeSet(PROP_DISABLE_CUTADAPT, false)==false;
			}
		
		private Integer getLane() {
			for(int i=0;i< this.fastqs.size();++i)
				{
				String fname= this.fastqs.get(i).getFilename();
				if(!fname.endsWith(".fastq.gz")) continue;
				String tokens[]=fname.split("[_]");
				for(int j=0;j+1< tokens.length;++j)
					{
					if(tokens[j].matches("L[0-9]+") && tokens[j+1].equals("R"+(i+1)))
						{
						return Integer.parseInt(tokens[j].substring(1).trim());
						}
					}
				
				}
			return null;
			}
		
		public  String bwamem()
			{
			if(getSample().isBamAlreadyProvided()) return "";
			final String inputs[]=new String[2];

			final StringBuilder sb= new StringBuilder().
				append(this.getSortedFilename()).
				append(":").
				append(this.get(0).getFilename());
			if(this.fastqs.size()>1)
				{
				sb.append(" ").
				append(this.get(1).getFilename());
				inputs[0]=this.get(0).getFilename();
				inputs[1]=this.get(1).getFilename();
				}
			
			if(!getAttribute(PROP_MAPPING_REGION,"").isEmpty())
				{
				sb.append(" ").append(getAttribute(PROP_MAPPING_REGION,""));
				}
			
			sb.append("\n").
				append(rulePrefix());
			
			if(this.fastqs.size()==1 && this.fastqs.get(0).isBam())
				{
				inputs[0]="$(addsuffix .bam2fq.R1.fq.gz,$@)";
				inputs[1]="$(addsuffix .bam2fq.R2.fq.gz,$@)";

				sb.append(" && $(java.exe)  -Djava.io.tmpdir=$(dir $@) -jar \"$(picard.jar)\" SamToFastq ").
					append(" I=$<  FASTQ=").append(inputs[0]).
					append(" SECOND_END_FASTQ=").append(inputs[1]).
					append(" UNPAIRED_FASTQ=$(addsuffix .bam2fq.unpaired.fq.gz,$@) ").
					append(" VALIDATION_STRINGENCY=SILENT VERBOSITY=INFO TMP_DIR=$(dir $@) ").
					append(" && rm --verbose $(addsuffix .bam2fq.unpaired.fq.gz,$@) ")
					;
				}
			
			String fq1;
			String fq2;
			if(isUsingCutadapt())
				{
				fq1="$(addsuffix .tmp.R1.fq,$@)";
				fq2="$(addsuffix .tmp.R2.fq,$@)";
				for(int side=0;side<2;++side)
					{
					sb.append(" && ${cutadapt.exe} -a ").
					append(getAttribute(side==0?PROP_CUTADAPT_ADAPTER5:PROP_CUTADAPT_ADAPTER3)).
					append(" ").
					append(inputs[side]).
					append(" -o $(addsuffix .tmp.fq,$@) > /dev/null ").
					append("&& awk '{if(NR%4==2 && length($$0)==0) { printf(\"N\\n\");} else if(NR%4==0 && length($$0)==0) { printf(\"#\\n\");} else {print;}}' ").
					append("$(addsuffix .tmp.fq,$@) > ").
					append(side==0?fq1:fq2).
					append(" && rm   --verbose $(addsuffix .tmp.fq,$@) ");
					}
				}
			else
				{
				fq1=inputs[0];
				fq2=inputs[1];
				}
			final Integer laneIndex= getLane();
			
			
			final String bwaRef="$(dir $(REF))index-bwa-0.7.12/$(notdir $(REF))";

			
			sb.append(" && $(bwa.exe) mem -t ").
				append(getAttribute(PROP_BWA_MEM_NTHREADS)).
				append(" -M -H '@CO\\tDate:"
					+ new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
				if(getProject().hasCapture())
					{	
					sb.append(" Capture was ").
						append(getProject().getCapture().getName()).
						append(" : ").
						append(getProject().getCapture().getBedFilename()).
						append(" ")
						;
					}
				sb.append(" Alignment of $(notdir $^) for project ").
					append(getProject().getName()+":"+getProject().getDescription()+"'").
				append(" -R '@RG\\tID:"+ getSample().getName()).
				append(laneIndex==null?"":"_L"+laneIndex).
				append("\\tLB:"+getSample()+"\\tSM:"+getSample() +"\\tPL:illumina\\tCN:Nantes' \""+bwaRef+"\" "+fq1+" "+fq2+"  |");
			
			if(!getAttribute(PROP_MAPPING_REGION,"").isEmpty())
				{
				final String rgn=getAttribute(PROP_MAPPING_REGION,"");
				sb.append("$(samtools.exe) view -b -u -T ${REF} -L \"").append(rgn).append("\" - | ");
				}

				
				
			sb.append("$(samtools.exe) sort --reference \"$(REF)\" -l ").
				append(getAttribute(PROP_DEFAULT_COMPRESSION_LEVEL)).
				append(" -@ ").
				append(getAttribute(PROP_SAMTOOLS_SORT_NTHREADS)).
				append(" -O bam -o $(addsuffix .tmp.bam,$@) -T  $(dir $@)"+getTmpPrefixToken()+getSample().getName()+".tmp_sort_"+getIndex()+" - && ");
			
			
			
			sb.append("mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\"")
				;
			if(isUsingCutadapt())
				{
				sb.append(" && rm --verbose -f "+fq1+" "+fq2);
				}
			if(this.fastqs.size()==1 && this.fastqs.get(0).isBam())
				{
				sb.append(" && rm --verbose -f  "+inputs[0]+" "+inputs[1]);
				}
			
			return sb.append("\n").
				toString();
			}
		
		@Override
		public String toString() {
			return getSample().getName()+".fastqs["+this.getIndex()+"]";
			}
		}
	
	private class Fastq extends HasAttributes
		{
		private final PairedFastq pair;
		private final String filename;
		Fastq( final PairedFastq pair,final JsonElement e ) throws IOException {
			super(e);
			this.pair=pair;
			filename= e.getAsString();
			}
		public String getFilename() {
			return filename;
			}
		public boolean isBam() {
			return getFilename().endsWith(".bam");
			}
		public boolean isFastq() {
			final String s=getFilename().toLowerCase();
			return s.endsWith(".fq") || s.endsWith(".fq.gz") || s.endsWith(".fastq") || s.endsWith(".fastq.gz");
			}
		@Override
		public String toString() {
			return getFilename();
			}
	
		public PairedFastq getPair() {
			return pair;
			}
		
		@Override
		public PairedFastq getParent() {
			return getPair();
			}
		}
	
	
	private PrintWriter out = new PrintWriter(System.out);
	
	public void execute(final Project project)
		{
		out.println("include ../../config/config.mk");
		out.println("OUT="+project.getOutputDirectory());
		
		
		out.println("# 1 : bam name");
		out.println("define bam_and_picardbai");
		out.println("$(1) $(call picardbai,$(1))");
		out.println("endef");
		
		out.println("#1 : bam filename");
		out.println("define picardbai");
		out.println("$(patsubst %.bam,%.bai,$(1))");
		out.println("endef");

		
		out.println("gatk.bundle.dir=$(dir $(REF))");
		out.println("gatk.bundle.dbsnp.vcf=$(gatk.bundle.dir)dbsnp_138.b37.vcf");
		out.println("gatk.bundle.onek.indels.vcf=$(gatk.bundle.dir)1000G_phase1.indels.b37.vcf");
		out.println("gatk.bundle.onek.snp.vcf=$(gatk.bundle.dir)1000G_phase1.snps.high_confidence.b37.vcf");
		out.println("gatk.bundle.mills.indels.vcf=$(gatk.bundle.dir)Mills_and_1000G_gold_standard.indels.b37.vcf");
		out.println("gatk.bundle.hapmap.vcf=$(gatk.bundle.dir)hapmap_3.3.b37.vcf");
		out.println("exac.vcf?=/commun/data/pubdb/broadinstitute.org/exac/0.3/ExAC.r0.3.sites.vep.vcf.gz");
		
		
		//out.println("define run_jvarkit\n${java.exe} -jar ${jvarkit.dir}/$(1).jar\nendef");

		if(!reannotatevcfs) 
			{
			out.print("all:"+project.getHapCallerAnnotationVcf());
			
			if(project.isTrue(PROP_USE_LUMPY_EXPRESS))
				{
				out.print(" ");
				out.print(project.getLumpyVcf());
				}
			
			out.println();
			
			out.println("all_final_bam:"+project.getSamples().stream().
					map(S->S.getFinalBamBai()).
					collect(Collectors.joining(" \\\n\t"))
					);
			
			out.println("all_recal_bam:"+project.getSamples().stream().
					filter(S->!S.isBamAlreadyProvided()).
					map(S->S.getRecalBamBai()).
					collect(Collectors.joining(" \\\n\t"))
					);
			
			out.println("all_realigned_bam:"+project.getSamples().stream().
					filter(S->!S.isBamAlreadyProvided()).
					map(S->S.getRealignBamBai()).
					collect(Collectors.joining(" \\\n\t"))
					);
			
			out.println("all_markdup_bam:"+project.getSamples().stream().
					filter(S->!S.isBamAlreadyProvided()).
					map(S->S.getMarkdupBamBai()).
					collect(Collectors.joining(" \\\n\t"))
					);
			
			out.println("all_merged_bam:"+project.getSamples().stream().
					filter(S->!S.isBamAlreadyProvided()).
					map(S->S.getMergedBamBai()).
					collect(Collectors.joining(" \\\n\t"))
					);
			
			
			project.lumpyExpress(out);
				
			
			
			for(final Sample sample: project.getSamples())
				{
				if(sample.getPairs().size()==1)
					{	
					out.println(sample.getPairs().get(0).bwamem());
					
					}
				else if(sample.getPairs().size()>1)
					{
					sample.getPairs().forEach(P->{out.println(P.bwamem());});
					out.println(sample.mergeSortedBams());
					}
				
				out.println(sample.markDuplicates());
				out.println(sample.realignAroundIndels());
				out.println(sample.recalibrateBam());
				out.println(sample.finalBam());
				}
			if(project.hasCapture())
				{
				project.getCapture().prepareCapture(out);
				}
			project.getPedigree().build(out);
			project.bamList(out);
			project.callVariants(out);
			project.annotateVariants(out);
			} 
		else
			/* reannotate vcf */ 
			{
			
			out.println("all:"+ project.vcfToReannotate.stream().map(S->project.getReannoteVcfFilename(S)).collect(Collectors.joining(" ")));
			project.reannoteVcfs(out);
			}
		
		out.println();
		out.flush();
		}
	
	private static JsonElement readJsonFile(final File jsonFile) throws IOException {
		IOUtil.assertFileIsReadable(jsonFile);
		FileReader r=new FileReader(jsonFile);
		JsonParser jsparser=new JsonParser();
		JsonElement root=jsparser.parse(r);
		r.close();
		return root;
		}
	@Override
	public int doWork(final List<String> args) {
		if(this.dumpAttributes)
			{
			for(final PropertyKey prop:NgsWorkflow.name2propertyKey.values())
				{
				System.out.print("\""+prop.getKey()+"\"\t"+prop.getDescription());
				if(prop.getDefaultValue()!=null && !prop.getDefaultValue().isEmpty())
					{
					System.out.print("\t\""+prop.getDefaultValue()+"\"");
					}
				System.out.println();
				}
			return 0;
			}
		try
			{
			if(args.size()!=1) throw new JvarkitException.CommandLineError("exepcted one and only one json file as input");
			final File jsonFile=new File(args.get(0));
			final JsonElement root=readJsonFile(jsonFile);
			final Project proj=new Project(root);
			execute(proj);
			out.flush();
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
	
	/* curl "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz" |gunzip -c | cut -f2,3,4 | grep -v _ | grep -v GL | grep -v MT | sed 's/^chr//' | sort -t '   ' -k1,1V -k2,2n -k3,3n | /commun/data/packages/bedtools/bedtools2-2.25.0/bin/complementBed -i - -g <(cut -f 1,2 /commun/data/pubdb/broadinstitute.org/bundle/1.5/b37/human_g1k_v37.fasta.fai) | sort -t '        ' -k1,1V -k2,2n -k3,3n | grep -v GL | awk -F '   ' '($2!=$3)' | awk '{printf("L.add(new IntervalSplit(%s,%s,%s));\n",$1,int($2)-1,int($3)+1);}' */
	public List<IntervalSplit> getIntervalSplitListGapForGrch37() {
		final List<IntervalSplit> L=new ArrayList<>();
		L.add(new IntervalSplit(1,10000,177417));
		L.add(new IntervalSplit(1,227416,267720));
		L.add(new IntervalSplit(1,317718,471369));
		L.add(new IntervalSplit(1,521367,2634221));
		L.add(new IntervalSplit(1,2684219,3845269));
		L.add(new IntervalSplit(1,3995267,13052999));
		L.add(new IntervalSplit(1,13102997,13219913));
		L.add(new IntervalSplit(1,13319911,13557163));
		L.add(new IntervalSplit(1,13607161,17125659));
		L.add(new IntervalSplit(1,17175657,29878083));
		L.add(new IntervalSplit(1,30028081,103863907));
		L.add(new IntervalSplit(1,103913905,120697157));
		L.add(new IntervalSplit(1,120747155,120936696));
		L.add(new IntervalSplit(1,121086694,121485435));
		L.add(new IntervalSplit(1,142535433,142731023));
		L.add(new IntervalSplit(1,142781021,142967762));
		L.add(new IntervalSplit(1,143117760,143292817));
		L.add(new IntervalSplit(1,143342815,143544526));
		L.add(new IntervalSplit(1,143644524,143771003));
		L.add(new IntervalSplit(1,143871001,144095784));
		L.add(new IntervalSplit(1,144145782,144224482));
		L.add(new IntervalSplit(1,144274480,144401745));
		L.add(new IntervalSplit(1,144451743,144622414));
		L.add(new IntervalSplit(1,144672412,144710725));
		L.add(new IntervalSplit(1,144810723,145833119));
		L.add(new IntervalSplit(1,145883117,146164651));
		L.add(new IntervalSplit(1,146214649,146253300));
		L.add(new IntervalSplit(1,146303298,148026039));
		L.add(new IntervalSplit(1,148176037,148361359));
		L.add(new IntervalSplit(1,148511357,148684148));
		L.add(new IntervalSplit(1,148734146,148954461));
		L.add(new IntervalSplit(1,149004459,149459646));
		L.add(new IntervalSplit(1,149509644,205922708));
		L.add(new IntervalSplit(1,206072706,206332222));
		L.add(new IntervalSplit(1,206482220,223747847));
		L.add(new IntervalSplit(1,223797845,235192212));
		L.add(new IntervalSplit(1,235242210,248908211));
		L.add(new IntervalSplit(1,249058209,249240622));
		L.add(new IntervalSplit(2,10000,3529312));
		L.add(new IntervalSplit(2,3579311,5018789));
		L.add(new IntervalSplit(2,5118787,16279725));
		L.add(new IntervalSplit(2,16329723,21153114));
		L.add(new IntervalSplit(2,21178112,87668207));
		L.add(new IntervalSplit(2,87718205,89630437));
		L.add(new IntervalSplit(2,89830435,90321526));
		L.add(new IntervalSplit(2,90371524,90545104));
		L.add(new IntervalSplit(2,91595102,92326172));
		L.add(new IntervalSplit(2,95326170,110109338));
		L.add(new IntervalSplit(2,110251336,149690583));
		L.add(new IntervalSplit(2,149790581,234003742));
		L.add(new IntervalSplit(2,234053740,239801979));
		L.add(new IntervalSplit(2,239831977,240784133));
		L.add(new IntervalSplit(2,240809131,243102477));
		L.add(new IntervalSplit(2,243152475,243189374));
		L.add(new IntervalSplit(3,60000,66170270));
		L.add(new IntervalSplit(3,66270269,90504855));
		L.add(new IntervalSplit(3,93504853,194041962));
		L.add(new IntervalSplit(3,194047250,197962431));
		L.add(new IntervalSplit(4,10000,1423146));
		L.add(new IntervalSplit(4,1478645,8799204));
		L.add(new IntervalSplit(4,8818202,9274643));
		L.add(new IntervalSplit(4,9324641,31820918));
		L.add(new IntervalSplit(4,31837416,32834639));
		L.add(new IntervalSplit(4,32840637,40296397));
		L.add(new IntervalSplit(4,40297095,49338942));
		L.add(new IntervalSplit(4,49488940,49660118));
		L.add(new IntervalSplit(4,52660116,59739334));
		L.add(new IntervalSplit(4,59789332,75427380));
		L.add(new IntervalSplit(4,75452278,191044277));
		L.add(new IntervalSplit(5,10000,17530657));
		L.add(new IntervalSplit(5,17580656,46405642));
		L.add(new IntervalSplit(5,49405640,91636129));
		L.add(new IntervalSplit(5,91686127,138787074));
		L.add(new IntervalSplit(5,138837072,155138728));
		L.add(new IntervalSplit(5,155188726,180905261));
		L.add(new IntervalSplit(6,60000,58087659));
		L.add(new IntervalSplit(6,58137658,58780167));
		L.add(new IntervalSplit(6,61880165,62128590));
		L.add(new IntervalSplit(6,62178588,95680544));
		L.add(new IntervalSplit(6,95830542,157559468));
		L.add(new IntervalSplit(6,157609466,157641301));
		L.add(new IntervalSplit(6,157691299,167942074));
		L.add(new IntervalSplit(6,168042072,170279973));
		L.add(new IntervalSplit(6,170329971,171055068));
		L.add(new IntervalSplit(7,10000,232484));
		L.add(new IntervalSplit(7,282483,50370632));
		L.add(new IntervalSplit(7,50410630,58054332));
		L.add(new IntervalSplit(7,61054330,61310514));
		L.add(new IntervalSplit(7,61360512,61460466));
		L.add(new IntervalSplit(7,61510464,61677021));
		L.add(new IntervalSplit(7,61727019,61917158));
		L.add(new IntervalSplit(7,61967156,74715725));
		L.add(new IntervalSplit(7,74765723,100556044));
		L.add(new IntervalSplit(7,100606042,130154524));
		L.add(new IntervalSplit(7,130254522,139379378));
		L.add(new IntervalSplit(7,139404376,142048196));
		L.add(new IntervalSplit(7,142098194,142276198));
		L.add(new IntervalSplit(7,142326196,143347898));
		L.add(new IntervalSplit(7,143397896,154270635));
		L.add(new IntervalSplit(7,154370633,159128664));
		L.add(new IntervalSplit(8,10000,7474649));
		L.add(new IntervalSplit(8,7524648,12091855));
		L.add(new IntervalSplit(8,12141853,43838888));
		L.add(new IntervalSplit(8,46838886,48130500));
		L.add(new IntervalSplit(8,48135598,86576452));
		L.add(new IntervalSplit(8,86726450,142766516));
		L.add(new IntervalSplit(8,142816514,145332589));
		L.add(new IntervalSplit(8,145432587,146304023));
		L.add(new IntervalSplit(9,10000,39663686));
		L.add(new IntervalSplit(9,39713685,39974797));
		L.add(new IntervalSplit(9,40024795,40233030));
		L.add(new IntervalSplit(9,40283028,40425835));
		L.add(new IntervalSplit(9,40475833,40940342));
		L.add(new IntervalSplit(9,40990340,41143215));
		L.add(new IntervalSplit(9,41193213,41365794));
		L.add(new IntervalSplit(9,41415792,42613956));
		L.add(new IntervalSplit(9,42663954,43213699));
		L.add(new IntervalSplit(9,43313697,43946570));
		L.add(new IntervalSplit(9,43996568,44676647));
		L.add(new IntervalSplit(9,44726645,44908294));
		L.add(new IntervalSplit(9,44958292,45250204));
		L.add(new IntervalSplit(9,45350202,45815522));
		L.add(new IntervalSplit(9,45865520,46216431));
		L.add(new IntervalSplit(9,46266429,46461040));
		L.add(new IntervalSplit(9,46561038,47060134));
		L.add(new IntervalSplit(9,47160132,47317680));
		L.add(new IntervalSplit(9,65467678,65918361));
		L.add(new IntervalSplit(9,65968359,66192216));
		L.add(new IntervalSplit(9,66242214,66404657));
		L.add(new IntervalSplit(9,66454655,66614196));
		L.add(new IntervalSplit(9,66664194,66863344));
		L.add(new IntervalSplit(9,66913342,67107835));
		L.add(new IntervalSplit(9,67207833,67366297));
		L.add(new IntervalSplit(9,67516295,67987999));
		L.add(new IntervalSplit(9,68137997,68514182));
		L.add(new IntervalSplit(9,68664180,68838947));
		L.add(new IntervalSplit(9,68988945,69278386));
		L.add(new IntervalSplit(9,69328384,70010543));
		L.add(new IntervalSplit(9,70060541,70218730));
		L.add(new IntervalSplit(9,70318728,70506536));
		L.add(new IntervalSplit(9,70556534,70735469));
		L.add(new IntervalSplit(9,70835467,92343417));
		L.add(new IntervalSplit(9,92443415,92528797));
		L.add(new IntervalSplit(9,92678795,133073061));
		L.add(new IntervalSplit(9,133223059,137041194));
		L.add(new IntervalSplit(9,137091192,139166998));
		L.add(new IntervalSplit(9,139216996,141153432));
		L.add(new IntervalSplit(10,60000,17974675));
		L.add(new IntervalSplit(10,18024674,38818836));
		L.add(new IntervalSplit(10,38868834,39154936));
		L.add(new IntervalSplit(10,42354934,42546688));
		L.add(new IntervalSplit(10,42596686,46426965));
		L.add(new IntervalSplit(10,46476963,47429170));
		L.add(new IntervalSplit(10,47529168,47792477));
		L.add(new IntervalSplit(10,47892475,48055708));
		L.add(new IntervalSplit(10,48105706,49095537));
		L.add(new IntervalSplit(10,49195535,51137411));
		L.add(new IntervalSplit(10,51187409,51398846));
		L.add(new IntervalSplit(10,51448844,125869473));
		L.add(new IntervalSplit(10,125919471,128616070));
		L.add(new IntervalSplit(10,128766068,133381405));
		L.add(new IntervalSplit(10,133431403,133677528));
		L.add(new IntervalSplit(10,133727526,135524748));
		L.add(new IntervalSplit(11,60000,1162759));
		L.add(new IntervalSplit(11,1212758,50783854));
		L.add(new IntervalSplit(11,51090852,51594206));
		L.add(new IntervalSplit(11,54694204,69089802));
		L.add(new IntervalSplit(11,69139800,69724696));
		L.add(new IntervalSplit(11,69774694,87688379));
		L.add(new IntervalSplit(11,87738377,96287585));
		L.add(new IntervalSplit(11,96437583,134946517));
		L.add(new IntervalSplit(12,60000,95739));
		L.add(new IntervalSplit(12,145738,7189877));
		L.add(new IntervalSplit(12,7239875,34856695));
		L.add(new IntervalSplit(12,37856693,109373471));
		L.add(new IntervalSplit(12,109423469,122530624));
		L.add(new IntervalSplit(12,122580622,132706993));
		L.add(new IntervalSplit(12,132806991,133841896));
		L.add(new IntervalSplit(13,19020000,86760324));
		L.add(new IntervalSplit(13,86910323,112353995));
		L.add(new IntervalSplit(13,112503993,114325994));
		L.add(new IntervalSplit(13,114425992,114639949));
		L.add(new IntervalSplit(13,114739947,115109879));
		L.add(new IntervalSplit(14,19000000,107289540));
		L.add(new IntervalSplit(15,20000000,20894633));
		L.add(new IntervalSplit(15,20935074,21398820));
		L.add(new IntervalSplit(15,21884999,22212115));
		L.add(new IntervalSplit(15,22262113,22596194));
		L.add(new IntervalSplit(15,22646192,23514854));
		L.add(new IntervalSplit(15,23564852,29159444));
		L.add(new IntervalSplit(15,29209442,82829646));
		L.add(new IntervalSplit(15,82879644,84984474));
		L.add(new IntervalSplit(15,85034472,102521393));
		L.add(new IntervalSplit(16,60000,8636921));
		L.add(new IntervalSplit(16,8686920,34023151));
		L.add(new IntervalSplit(16,34173149,35285802));
		L.add(new IntervalSplit(16,46385800,88389384));
		L.add(new IntervalSplit(16,88439382,90294754));
		L.add(new IntervalSplit(17,-1,296627));
		L.add(new IntervalSplit(17,396625,21566609));
		L.add(new IntervalSplit(17,21666607,22263007));
		L.add(new IntervalSplit(17,25263005,34675849));
		L.add(new IntervalSplit(17,34725847,62410761));
		L.add(new IntervalSplit(17,62460759,77546462));
		L.add(new IntervalSplit(17,77596460,79709050));
		L.add(new IntervalSplit(17,79759048,81195211));
		L.add(new IntervalSplit(18,10000,15410898));
		L.add(new IntervalSplit(18,18510897,52059137));
		L.add(new IntervalSplit(18,52209135,72283354));
		L.add(new IntervalSplit(18,72333352,75721821));
		L.add(new IntervalSplit(18,75771819,78017249));
		L.add(new IntervalSplit(19,60000,7346004));
		L.add(new IntervalSplit(19,7396003,8687199));
		L.add(new IntervalSplit(19,8737197,20523416));
		L.add(new IntervalSplit(19,20573414,24631783));
		L.add(new IntervalSplit(19,27731781,59118984));
		L.add(new IntervalSplit(20,60000,26319569));
		L.add(new IntervalSplit(20,29419568,29653909));
		L.add(new IntervalSplit(20,29803907,34897086));
		L.add(new IntervalSplit(20,34947084,61091438));
		L.add(new IntervalSplit(20,61141436,61213370));
		L.add(new IntervalSplit(20,61263368,62965521));
		L.add(new IntervalSplit(21,9411193,9595548));
		L.add(new IntervalSplit(21,9645547,9775438));
		L.add(new IntervalSplit(21,9825436,10034921));
		L.add(new IntervalSplit(21,10084919,10215977));
		L.add(new IntervalSplit(21,10365975,10647897));
		L.add(new IntervalSplit(21,10697895,11188130));
		L.add(new IntervalSplit(21,14338128,42955560));
		L.add(new IntervalSplit(21,43005558,44632665));
		L.add(new IntervalSplit(21,44682663,48119896));
		L.add(new IntervalSplit(22,16050000,16697850));
		L.add(new IntervalSplit(22,16847849,20509432));
		L.add(new IntervalSplit(22,20609430,50364778));
		L.add(new IntervalSplit(22,50414776,51244567));

		
		L.add(new IntervalSplit("X",60000,94821));
		L.add(new IntervalSplit("X",144820,231385));
		L.add(new IntervalSplit("X",281383,1047558));
		L.add(new IntervalSplit("X",1097556,1134114));
		L.add(new IntervalSplit("X",1184112,1264235));
		L.add(new IntervalSplit("X",1314233,2068239));
		L.add(new IntervalSplit("X",2118237,7623883));
		L.add(new IntervalSplit("X",7673881,10738675));
		L.add(new IntervalSplit("X",10788673,37098257));
		L.add(new IntervalSplit("X",37148255,49242998));
		L.add(new IntervalSplit("X",49292996,49974174));
		L.add(new IntervalSplit("X",50024172,52395915));
		L.add(new IntervalSplit("X",52445913,58582013));
		L.add(new IntervalSplit("X",61682011,76653693));
		L.add(new IntervalSplit("X",76703691,113517669));
		L.add(new IntervalSplit("X",113567667,115682291));
		L.add(new IntervalSplit("X",115732289,120013236));
		L.add(new IntervalSplit("X",120063234,143507325));
		L.add(new IntervalSplit("X",143557323,148906425));
		L.add(new IntervalSplit("X",148956423,149032063));
		L.add(new IntervalSplit("X",149082061,152277100));
		L.add(new IntervalSplit("X",152327098,155260561));
		L.add(new IntervalSplit("Y",10000,44821));
		L.add(new IntervalSplit("Y",94820,181385));
		L.add(new IntervalSplit("Y",231383,997558));
		L.add(new IntervalSplit("Y",1047556,1084114));
		L.add(new IntervalSplit("Y",1134112,1214235));
		L.add(new IntervalSplit("Y",1264233,2018239));
		L.add(new IntervalSplit("Y",2068237,8914956));
		L.add(new IntervalSplit("Y",8964954,9241323));
		L.add(new IntervalSplit("Y",9291321,10104554));
		L.add(new IntervalSplit("Y",13104552,13143955));
		L.add(new IntervalSplit("Y",13193953,13748579));
		L.add(new IntervalSplit("Y",13798577,20143886));
		L.add(new IntervalSplit("Y",20193884,22369680));
		L.add(new IntervalSplit("Y",22419678,23901429));
		L.add(new IntervalSplit("Y",23951427,28819362));
		L.add(new IntervalSplit("Y",58819360,58917657));
		L.add(new IntervalSplit("Y",58967655,59363567));
		return L;
		}
	
	public static void main(String[] args) {
		new NgsWorkflow().instanceMainWithExit(args);
	}
	
}
