package com.github.lindenb.jvarkit.tools.workflow;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

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
			def("1").
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
		protected String indexBam(final String indexFile,final String bamFile) {
			return new StringBuilder(indexFile).append(bamFile).append("\n\t${samtools.exe} index $<").toString();
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
		
		public List<Sample> getSamples() { return this._samples;}
		@Override  public HasAttributes getParent() { return null;}
		
		public  String getOutputDirectory() {
			return getAttribute(PROP_PROJECT_OUTPUT_DIRECTORY);
			}
		public  String getSamplesDirectory() {
			return "${OUT}/Samples";
			}
	
		}
	
	
	private class Sample extends HasAttributes
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
					final PairedFastq pair =new PairedFastq(this,this.pairs.size(),pairjson);
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
					append("$(java.exe) -Xmx3g -Djava.io.tmpdir=$(dir $@) -jar \"$(picard.jar)\" MarkDuplicates I=$< O=$(addsuffix .tmp.bam,$@) M=$(addsuffix .metrics,$@) ").
					append(" SO=coordinate AS=true CREATE_INDEX=false  COMPRESSION_LEVEL=9 VALIDATION_STRINGENCY=SILENT VERBOSITY=INFO TMP_DIR=$(dir $@)  ").
					append(" && mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\" ").
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
				return getDirectory()+"/${TMPPREFIX}."+getSample().getName()+".sorted.bam";
				}
			}
		private boolean isUsingCutadapt()
			{
			return true;
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
			sb.append("\n").
				append("\tmkdir -p $(dir $@  ) ");
			
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
			
			sb.append(" && $(bwa.exe) mem -t ").
				append(getAttribute(PROP_BWA_MEM_NTHREADS)).
				append(" -M -H '@CO\\tDate:"
					+ new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date())
					+" Alignment of $(notdir $^) for project "+
					getProject().getName()+":"+getProject().getDescription()+"'").
				append(" -R '@RG\\tID:"+ getSample().getName()).
				append(laneIndex==null?"":"_L"+laneIndex).
				append("\\tLB:"+getSample()+"\\tSM:"+getSample() +"\\tPL:illumina\\tCN:Nantes' ${REF} "+fq1+" "+fq2+"  |").
				append("$(samtools.exe) sort --reference $(REF) -l 9 -@ ").
				append(getAttribute(PROP_SAMTOOLS_SORT_NTHREADS)).
				append(" -O bam -o $(addsuffix .tmp.bam,$@) -T  $(dir $@)${TMPPREFIX}_sort_"+getIndex()+" - && ").
				append("mv --verbose \"$(addsuffix .tmp.bam,$@)\" \"$@\"")
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
		}
	
	private class Fastq extends HasAttributes
		{
		final PairedFastq pair;
		final String filename;
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
		out.println("include ../../config/config.mk");
		out.println("OUT="+project.getOutputDirectory());
		
		out.println("all_merged_bam:"+project.getSamples().stream().
					map(S->S.getMergedBam()).
					collect(Collectors.joining(" \\\n\t"))
					);
		
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
			Project proj=new Project(root);
			execute(proj);
			out.flush();
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
