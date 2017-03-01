package com.github.lindenb.jvarkit.tools.workflow;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

public class NgsWorkflow extends AbstractNgsWorkflow
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(NgsWorkflow.class);
	
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
	
	
	private static final PropertyKey PROP_PROJECT_OUTPUT_DIRECTORY = key("project.output.directory").
			description("Project output directory").
			def("${OUTDIR}").
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
		
		public boolean isAttributeSet(final PropertyKey key,boolean def)
			{
			final String s= getAttribute(key,def?"true":"false");
			if(s==null) throw new RuntimeException("key "+key+" is not defined for "+this.toString());
			if(s.equals("1") || s.equals("true") || s.equals("yes") || s.equals("T") || s.equals("Y")) return true;
			if(s.equals("0") || s.equals("false") || s.equals("no") || s.equals("F") || s.equals("N")) return false;
			throw new RuntimeException("bad boolean value for "+key+" : "+s);
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
		protected void indexBam(final String indexFile,final String bamFile) {
			out.println(indexFile+":"+bamFile);
			out.println("\t${samtools.exe} index $<");
			}
		}
	
	/** describe a NGS project */
	private class Project extends HasAttributes
		{
		private final String name;
		private final String description;
		private final List<Sample> _samples=new ArrayList<>();
		Project(final JsonElement root) throws IOException
			{
			super(root);
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
			
			final Set<String> seen=new HashSet<>(); 
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
		
		public String getName() {
			return name;
			}
		
		public String getDescription() {
			return description;
		}
		
		public List<Sample> getSamples() { return this._samples;}
		@Override  public HasAttributes getParent() { return null;}
		
		public  String getOutputDirectory() {
			return getAttribute(PROP_PROJECT_OUTPUT_DIRECTORY);
			}
		public  String getSamplesDirectory() {
			return getOutputDirectory()+"/Samples";
			}
	
		}
	
	
	class Sample extends HasAttributes
		{
		private final Project project;
		private final String name;
		private final List<PairedFastq> pairs;
		Sample( final Project project,final JsonObject root) throws IOException
			{
			super(root);
			if(!root.isJsonObject()) throw new IOException("sample json is not object");
			final JsonObject json = root.getAsJsonObject();
	
			this.project = project;
			
			if(!json.has("name")) throw new IOException("@name missing in sample");
			this.name= json.get("name").getAsString();
			if(json.has("fastqs"))
				{
				this.pairs = new ArrayList<>();
				for(final JsonElement pairjson :json.get("fastqs").getAsJsonArray())
					{
					PairedFastq pair =new PairedFastq(this,this.pairs.size(),pairjson.getAsJsonObject());
					this.pairs.add(pair);
					}
				}
			else
				{
				this.pairs = Collections.emptyList();
				}
			}
		
		public Project getProject() { return this.project;}
		public List<PairedFastq> getPairs() { return this.pairs;}
		public String getName() { return this.name;}
		
		@Override public Project getParent() { return getProject();}
		
		public String getDirectory() {
			return getProject().getSamplesDirectory()+"/"+getName();
			}
		
		public String getMergedBam() {
			return getDirectory()+"/"+getTmpPrefix()+getName()+".merged.bam";
			}
		public String getMergedBamBai() {
			return "$(addsuffix .bai,"+getMergedBam()+")";
			}
		public String getMarkdupBam() {
			return getDirectory()+"/"+getTmpPrefix()+getName()+".markdup.bam";
			}
	
		
		public String mergeSortedBams()
			{
			return new StringBuilder().
				append(getMergedBam()+" : "+getPairs().stream().map(P->P.getSortedFilename()).collect(Collectors.joining(" "))+"\n").
				append("\tmkdir -p $(dir $@) && ").
				append("$(java.exe)  -Djava.io.tmpdir=$(dir $@) -jar \"$(picard.jar)\" MergeSamFiles O=$(addsuffix .tmp.bam,$@) ").
				append(" SO=coordinate AS=true CREATE_INDEX=false  COMPRESSION_LEVEL=9 VALIDATION_STRINGENCY=SILENT USE_THREADING=true VERBOSITY=INFO TMP_DIR=$(dir $@) COMMENT=\"Merged from $(words $^) files.\" ").
				append("$(foreach B,$(filter %.bam,$^), I=$(B) ) && ").
				append("mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\" ").
				append("\n").
				toString();
			}
		
		public String markDuplicates()
			{
			return new StringBuilder().
					append(getMarkdupBam()+" : "+getMergedBam()+"\n").
					append("\tmkdir -p $(dir $@) && ").
					append("$(java.exe)  -Djava.io.tmpdir=$(dir $@) -jar \"$(picard.jar)\" MergeSamFiles O=$(addsuffix .tmp.bam,$@) ").
					append(" SO=coordinate AS=true CREATE_INDEX=false  COMPRESSION_LEVEL=9 VALIDATION_STRINGENCY=SILENT USE_THREADING=true VERBOSITY=INFO TMP_DIR=$(dir $@) COMMENT=\"Merged from $(words $^) files.\" ").
					append("$(foreach B,$(filter %.bam,$^), I=$(B) ) && ").
					append("mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\" ").
					append("\n").
					toString();
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
			if(root.isJsonObject())
				{
				final JsonObject json=root.getAsJsonObject();
				for(final JsonElement fqjson :json.get("fastqs").getAsJsonArray())
					{		
					Fastq fq =new Fastq(this,fqjson);
					this.fastqs.add(fq);
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
				return getDirectory()+"/${TMPPREFIX}."+getSample().getName()+".sorted.bam";
				}
			}
		
		public  String bwamem()
			{
			return new StringBuilder().
				append(this.getSortedFilename()).
				append(":").
				append(this.get(0)).
				append(" ").
				append(this.get(1)).
				append("\n").
				append("\tmkdir -p $(dir $@  ) && ").
				append(" $(bwa.exe) mem -t 1 -M -H '@CO\\tAlignment de $(notdir $^) par Pierre Lindenbaum pour Cedric LeCaignec 20170228' ").
				append(" -R '@RG\\tID:"+ getSample().getName()+"_"+getIndex()+"\\tLB:"+getSample()+"\\tSM:"+getSample() +"\\tPL:illumina\\tCN:Nantes' ${REF} $(word 1,$^)  $(word 2,$^)  |").
				append("$(samtools.exe) sort --reference $(REF) -l 9 -@ 1 -O bam -o $(addsuffix .tmp.bam,$@) -T  $(dir $@)${TMPPREFIX}_sort_"+getIndex()+" - && ").
				append("mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\"").
				append("\n").
				toString();
			}
		}
	
	class Fastq
		{
		final PairedFastq pair;
		final String filename;
		Fastq( final PairedFastq pair,final JsonElement e ) {
			this.pair=pair;
			filename= e.getAsString();
			}
		public String getFilename() {
			return filename;
			}
		}
	
	BiFunction<String, List<String>, String> mergeBamsFunctions = new BiFunction<String, List<String>, String>() {
		@Override
		public String apply(final String out, final List<String> list) {
			return new StringBuilder().
			append(out+" : "+String.join(" ", list)+"\n").
			append("\tmkdir -p $(dir $@) && ").
			append("$(java.exe)  -Djava.io.tmpdir=$(dir $@) -jar \"$(picard.jar)\" MergeSamFiles O=$(addsuffix .tmp.bam,$@) ").
			append(" SO=coordinate AS=true CREATE_INDEX=false  COMPRESSION_LEVEL=9 VALIDATION_STRINGENCY=SILENT USE_THREADING=true VERBOSITY=INFO TMP_DIR=$(dir $@) COMMENT=\"Merged from $(words $^) files.\" ").
			append("$(foreach B,$(filter %.bam,$^), I=$(B) ) && ").
			append("mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\" ").
			append("\n").
			toString();
		}
	};
	
	private PrintWriter out = new PrintWriter(System.out);
	
	public void execute(final Project project)
		{
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
			}
		
		}
	@Override
	public Collection<Throwable> call() throws Exception {
		if(super.dumpAttributes)
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
			return RETURN_OK;
			}
		FileReader r=null;
		try
			{
			final List<String> args=this.getInputFiles();
			if(args.size()!=1) return wrapException("exepcted one and only one json file as input");
			final File jsonFile=new File(args.get(0));
			IOUtil.assertFileIsReadable(jsonFile);
			r=new FileReader(jsonFile);
			JsonParser jsparser=new JsonParser();
			JsonElement root=jsparser.parse(r);
			r.close();
			Project proj=new Project(root.getAsJsonObject());
			execute(proj);
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(r);
			}
		}
	
	public static void main(String[] args) {
		new NgsWorkflow().instanceMainWithExit(args);
	}
	
}
