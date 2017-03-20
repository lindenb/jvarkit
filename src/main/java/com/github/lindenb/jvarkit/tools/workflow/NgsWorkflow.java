package com.github.lindenb.jvarkit.tools.workflow;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonPrimitive;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;

public class NgsWorkflow extends AbstractNgsWorkflow
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(NgsWorkflow.class);
	
	private enum RefSplitType {WHOLE_GENOME,WHOLE_CONTIG,INTERVAL};
	

	
	
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
		
		public  String getHapCallerGenotypedVcf() {
			return getVcfDirectory()+"/"+getTmpPrefix()+"Genotyped.vcf.gz";
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
		
		public StringBuilder haplotypeCaller() {
			final StringBuilder w=new StringBuilder();
			final List<String> vcfParts=new ArrayList<>();
			final List<RefSplit> refSplits = getHaplotypeCallerSplits();
			if(refSplits.isEmpty()) throw new IllegalArgumentException();
			
			for(final RefSplit split: refSplits)
				{
				final String vcfPart;
				switch(split.getType()) {
					case WHOLE_GENOME:  vcfPart= getHapCallerGenotypedVcf();break;
					default:  vcfPart= "$(addsuffix "+split.getToken()+".vcf.gz,"+getHapCallerGenotypedVcf()+")";
					}
				
				vcfParts.add(vcfPart);
				
				w.append(vcfPart).append(":").append(getFinalBamList());
				w.append(" ").append(getPedigree().getPedFilename());
				if( this.hasCapture())
	            	{
					w.append(" ").append(getCapture().getExtendedFilename());
	            	}
				w.append("\n");
				w.append(rulePrefix()+" && ");
				
				if( this.hasCapture())
	            	{
					switch(split.getType())
						{
						case WHOLE_GENOME: break;
						case INTERVAL:{
								IntervalSplit tmp= IntervalSplit.class.cast(split);
								w.append(" awk -F '\t' '($$1==\""+tmp.getInterval().getContig()+" && !("+ 
										tmp.getInterval().getEnd()+" < int($$2) || int($$3) <"+
										tmp.getInterval().getStart()+")' "+getCapture().getExtendedFilename()+" > $(addsuffix .bed,$@) && "); 
								break;
								}
						case WHOLE_CONTIG:
								{
								ContigSplit tmp= ContigSplit.class.cast(split);
								w.append(" awk -F '\t' '($$1==\""+tmp.getContig()+"\")' "+getCapture().getExtendedFilename()+" > $(addsuffix .bed,$@) && ");
								break;
								}
						default: throw new IllegalStateException();
						}
	            	}
				
				w.append(" $(java.exe)  -XX:ParallelGCThreads=5 -Xmx2g  -Djava.io.tmpdir=$(dir $@) -jar $(gatk.jar) -T HaplotypeCaller ");
				w.append("	-R $(REF) ");
				w.append("	--validation_strictness LENIENT ");
				w.append("	-I $< -o \"$(addsuffix .tmp.vcf.gz,$@)\" ");
				//w.print("	--num_cpu_threads_per_data_thread "+  this.getIntProperty("base-recalibrator-nct",1));
				w.append("	-l INFO ");
				w.append("	-nct ").append(getAttribute(PROP_GATK_HAPCALLER_NCT));
				
				w.append("	--dbsnp \"$(gatk.bundle.dbsnp.vcf)\" ");
				w.append(" $(foreach A, PossibleDeNovo AS_FisherStrand AlleleBalance AlleleBalanceBySample BaseCountsBySample GCContent ClippingRankSumTest , --annotation ${A} ) ");
				w.append(" --pedigree ").append(getPedigree().getPedFilename());
				if( this.hasCapture())
		            {
					switch(split.getType())
						{
						case WHOLE_GENOME: 	w.append(" -L:"+getCapture().getName()+",BED "+getCapture().getExtendedFilename());
						case INTERVAL: //through...
						case WHOLE_CONTIG: w.append(" -L $(addsuffix .bed,$@) ");break;
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
				w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz,$@)\" \"$@\" ");
				w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz.tbi,$@)\" \"$(addsuffix .tbi,$@)\" ");
				
				if( this.hasCapture())
	            	{
					w.append(" && rm --verbose \"$(addsuffix .bed,$@)\" ");
	            	}
				w.append("\n");
				}
			
			if(vcfParts.size()>1) {
				w.append(getHapCallerGenotypedVcf()).append(":");
				for(final String vcfPart:vcfParts) {
					w.append(" \\\n\t").append(vcfPart);
					}
				w.append("\n");
				
				w.append(rulePrefix()+" && ${java.exe}   -Djava.io.tmpdir=$(dir $@)  -jar ${gatk.jar}  -T CombineVariants -R $(REF) "
						+ " -o $(addsuffix .tmp.vcf.gz,$@) -genotypeMergeOptions UNSORTED "
						);
				for(int k=0;k<vcfParts.size();++k) {
					w.append(" --variant:v").append(k).append(" ").append(vcfParts.get(k)).append(" ");
					}

				w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz,$@)\" \"$@\" ");
				w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz.tbi,$@)\" \"$(addsuffix .tbi,$@)\" ");
				w.append("\n");
				}
			
			return w;
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
			L.add(new SamtoolsCaller());
			return L;
			}
		
		@Override
		public String toString() {
			return getName();
			}
		
		
		
		
		private abstract class AbstractCaller
			{
			Project getProject() { return Project.this;}
			abstract List<RefSplit> getCallSplits();
			abstract String getTargetVcfFilename();
			
			abstract void call(final StringBuilder w,RefSplit split);
			
			void  print(final StringBuilder w) {
				final List<String> vcfParts=new ArrayList<>();
				final List<RefSplit> refSplits = this.getCallSplits();
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
									IntervalSplit tmp= IntervalSplit.class.cast(split);
									w.append(" awk -F '\t' '($$1==\""+tmp.getInterval().getContig()+" && !("+ 
											tmp.getInterval().getEnd()+" < int($$2) || int($$3) <"+
											tmp.getInterval().getStart()+")' "+getCapture().getExtendedFilename()+" > $(addsuffix .bed,$@) && "); 
									break;
									}
							case WHOLE_CONTIG:
									{
									ContigSplit tmp= ContigSplit.class.cast(split);
									w.append(" awk -F '\t' '($$1==\""+tmp.getContig()+"\")' "+getCapture().getExtendedFilename()+" > $(addsuffix .bed,$@) && ");
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
					
					w.append(rulePrefix()+" && ${java.exe}   -Djava.io.tmpdir=$(dir $@)  -jar ${gatk.jar}  -T CombineVariants -R $(REF) "
							+ " -o $(addsuffix .tmp.vcf.gz,$@) -genotypeMergeOptions UNSORTED "
							);
					for(int k=0;k<vcfParts.size();++k) {
						w.append(" --variant:v").append(k).append(" ").append(vcfParts.get(k)).append(" ");
						}
	
					w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz,$@)\" \"$@\" ");
					w.append(" && mv --verbose \"$(addsuffix .tmp.vcf.gz.tbi,$@)\" \"$(addsuffix .tbi,$@)\" ");
					w.append("\n");
					}
				}
			}
	
	
	private class UnifiedGenotyperCaller extends AbstractCaller
		{	
		@Override
		List<RefSplit> getCallSplits() {
			return getProject().getHaplotypeCallerSplits();
			}
		@Override String getTargetVcfFilename() { return  getProject().getHapCallerGenotypedVcf();}
		
		@Override void call(final StringBuilder w,final RefSplit split)
			{
			w.append(" $(java.exe)  -XX:ParallelGCThreads=5 -Xmx2g  -Djava.io.tmpdir=$(dir $@) -jar $(gatk.jar) -T HaplotypeCaller ");
			w.append("	-R $(REF) ");
			w.append("	--validation_strictness LENIENT ");
			w.append("	-I $< -o \"$(addsuffix .tmp.vcf.gz,$@)\" ");
			//w.print("	--num_cpu_threads_per_data_thread "+  this.getIntProperty("base-recalibrator-nct",1));
			w.append("	-l INFO ");
			w.append("	-nct ").append(getAttribute(PROP_GATK_HAPCALLER_NCT));
			
			w.append("	--dbsnp \"$(gatk.bundle.dbsnp.vcf)\" ");
			w.append(" $(foreach A, PossibleDeNovo AS_FisherStrand AlleleBalance AlleleBalanceBySample BaseCountsBySample GCContent ClippingRankSumTest , --annotation ${A} ) ");
			w.append(" --pedigree ").append(getProject().getPedigree().getPedFilename());
			if( getProject().hasCapture())
	            {
				switch(split.getType())
					{
					case WHOLE_GENOME: 	w.append(" -L:"+getProject().getCapture().getName()+",BED "+getCapture().getExtendedFilename());
					case INTERVAL: //through...
					case WHOLE_CONTIG: w.append(" -L $(addsuffix .bed,$@) ");break;
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
		
		@Override void call(final StringBuilder w,final RefSplit split)
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
				append("\\tLB:"+getSample()+"\\tSM:"+getSample() +"\\tPL:illumina\\tCN:Nantes' \""+bwaRef+"\" "+fq1+" "+fq2+"  |").
				append("$(samtools.exe) sort --reference \"$(REF)\" -l ").
				append(getAttribute(PROP_DEFAULT_COMPRESSION_LEVEL)).
				append(" -@ ").
				append(getAttribute(PROP_SAMTOOLS_SORT_NTHREADS)).
				append(" -O bam -o $(addsuffix .tmp.bam,$@) -T  $(dir $@)"+getTmpPrefixToken()+getSample().getName()+".tmp_sort_"+getIndex()+" - && ").
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
		
		out.println("all:"+project.getHapCallerGenotypedVcf());
		
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
		out.println(project.haplotypeCaller());
		
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
		try
			{
			final List<String> args=this.getInputFiles();
			if(args.size()!=1) return wrapException("exepcted one and only one json file as input");
			final File jsonFile=new File(args.get(0));
			JsonElement root=readJsonFile(jsonFile);
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
			}
		}
	
	
	public List<IntervalSplit> getIntervalSplitListGapForGrch37() {
		final List<IntervalSplit> L=new ArrayList<>();
		L.add(new IntervalSplit(1,10001,177416));
		L.add(new IntervalSplit(1,227417,267719));
		L.add(new IntervalSplit(1,317719,471368));
		L.add(new IntervalSplit(1,521368,2634220));
		L.add(new IntervalSplit(1,2684220,3845268));
		L.add(new IntervalSplit(1,3995268,13052998));
		L.add(new IntervalSplit(1,13102998,13219912));
		L.add(new IntervalSplit(1,13319912,13557162));
		L.add(new IntervalSplit(1,13607162,17125658));
		L.add(new IntervalSplit(1,17175658,29878082));
		L.add(new IntervalSplit(1,30028082,103863906));
		L.add(new IntervalSplit(1,103913906,120697156));
		L.add(new IntervalSplit(1,120747156,120936695));
		L.add(new IntervalSplit(1,121086695,121485434));
		L.add(new IntervalSplit(1,142535434,142731022));
		L.add(new IntervalSplit(1,142781022,142967761));
		L.add(new IntervalSplit(1,143117761,143292816));
		L.add(new IntervalSplit(1,143342816,143544525));
		L.add(new IntervalSplit(1,143644525,143771002));
		L.add(new IntervalSplit(1,143871002,144095783));
		L.add(new IntervalSplit(1,144145783,144224481));
		L.add(new IntervalSplit(1,144274481,144401744));
		L.add(new IntervalSplit(1,144451744,144622413));
		L.add(new IntervalSplit(1,144672413,144710724));
		L.add(new IntervalSplit(1,144810724,145833118));
		L.add(new IntervalSplit(1,145883118,146164650));
		L.add(new IntervalSplit(1,146214650,146253299));
		L.add(new IntervalSplit(1,146303299,148026038));
		L.add(new IntervalSplit(1,148176038,148361358));
		L.add(new IntervalSplit(1,148511358,148684147));
		L.add(new IntervalSplit(1,148734147,148954460));
		L.add(new IntervalSplit(1,149004460,149459645));
		L.add(new IntervalSplit(1,149509645,205922707));
		L.add(new IntervalSplit(1,206072707,206332221));
		L.add(new IntervalSplit(1,206482221,223747846));
		L.add(new IntervalSplit(1,223797846,235192211));
		L.add(new IntervalSplit(1,235242211,248908210));
		L.add(new IntervalSplit(1,249058210,249240621));
		L.add(new IntervalSplit(2,10001,3529311));
		L.add(new IntervalSplit(2,3579312,5018788));
		L.add(new IntervalSplit(2,5118788,16279724));
		L.add(new IntervalSplit(2,16329724,21153113));
		L.add(new IntervalSplit(2,21178113,87668206));
		L.add(new IntervalSplit(2,87718206,89630436));
		L.add(new IntervalSplit(2,89830436,90321525));
		L.add(new IntervalSplit(2,90371525,90545103));
		L.add(new IntervalSplit(2,91595103,92326171));
		L.add(new IntervalSplit(2,95326171,110109337));
		L.add(new IntervalSplit(2,110251337,149690582));
		L.add(new IntervalSplit(2,149790582,234003741));
		L.add(new IntervalSplit(2,234053741,239801978));
		L.add(new IntervalSplit(2,239831978,240784132));
		L.add(new IntervalSplit(2,240809132,243102476));
		L.add(new IntervalSplit(2,243152476,243189373));
		L.add(new IntervalSplit(3,60001,66170269));
		L.add(new IntervalSplit(3,66270270,90504854));
		L.add(new IntervalSplit(3,93504854,194041961));
		L.add(new IntervalSplit(3,194047251,197962430));
		L.add(new IntervalSplit(4,10001,1423145));
		L.add(new IntervalSplit(4,1478646,8799203));
		L.add(new IntervalSplit(4,8818203,9274642));
		L.add(new IntervalSplit(4,9324642,31820917));
		L.add(new IntervalSplit(4,31837417,32834638));
		L.add(new IntervalSplit(4,32840638,40296396));
		L.add(new IntervalSplit(4,40297096,49338941));
		L.add(new IntervalSplit(4,49488941,49660117));
		L.add(new IntervalSplit(4,52660117,59739333));
		L.add(new IntervalSplit(4,59789333,75427379));
		L.add(new IntervalSplit(4,75452279,191044276));
		L.add(new IntervalSplit(5,10001,17530656));
		L.add(new IntervalSplit(5,17580657,46405641));
		L.add(new IntervalSplit(5,49405641,91636128));
		L.add(new IntervalSplit(5,91686128,138787073));
		L.add(new IntervalSplit(5,138837073,155138727));
		L.add(new IntervalSplit(5,155188727,180905260));
		L.add(new IntervalSplit(6,60001,58087658));
		L.add(new IntervalSplit(6,58137659,58780166));
		L.add(new IntervalSplit(6,61880166,62128589));
		L.add(new IntervalSplit(6,62178589,95680543));
		L.add(new IntervalSplit(6,95830543,157559467));
		L.add(new IntervalSplit(6,157609467,157641300));
		L.add(new IntervalSplit(6,157691300,167942073));
		L.add(new IntervalSplit(6,168042073,170279972));
		L.add(new IntervalSplit(6,170329972,171055067));
		L.add(new IntervalSplit(7,10001,232483));
		L.add(new IntervalSplit(7,282484,50370631));
		L.add(new IntervalSplit(7,50410631,58054331));
		L.add(new IntervalSplit(7,61054331,61310513));
		L.add(new IntervalSplit(7,61360513,61460465));
		L.add(new IntervalSplit(7,61510465,61677020));
		L.add(new IntervalSplit(7,61727020,61917157));
		L.add(new IntervalSplit(7,61967157,74715724));
		L.add(new IntervalSplit(7,74765724,100556043));
		L.add(new IntervalSplit(7,100606043,130154523));
		L.add(new IntervalSplit(7,130254523,139379377));
		L.add(new IntervalSplit(7,139404377,142048195));
		L.add(new IntervalSplit(7,142098195,142276197));
		L.add(new IntervalSplit(7,142326197,143347897));
		L.add(new IntervalSplit(7,143397897,154270634));
		L.add(new IntervalSplit(7,154370634,159128663));
		L.add(new IntervalSplit(8,10001,7474648));
		L.add(new IntervalSplit(8,7524649,12091854));
		L.add(new IntervalSplit(8,12141854,43838887));
		L.add(new IntervalSplit(8,46838887,48130499));
		L.add(new IntervalSplit(8,48135599,86576451));
		L.add(new IntervalSplit(8,86726451,142766515));
		L.add(new IntervalSplit(8,142816515,145332588));
		L.add(new IntervalSplit(8,145432588,146304022));
		L.add(new IntervalSplit(9,10001,39663685));
		L.add(new IntervalSplit(9,39713686,39974796));
		L.add(new IntervalSplit(9,40024796,40233029));
		L.add(new IntervalSplit(9,40283029,40425834));
		L.add(new IntervalSplit(9,40475834,40940341));
		L.add(new IntervalSplit(9,40990341,41143214));
		L.add(new IntervalSplit(9,41193214,41365793));
		L.add(new IntervalSplit(9,41415793,42613955));
		L.add(new IntervalSplit(9,42663955,43213698));
		L.add(new IntervalSplit(9,43313698,43946569));
		L.add(new IntervalSplit(9,43996569,44676646));
		L.add(new IntervalSplit(9,44726646,44908293));
		L.add(new IntervalSplit(9,44958293,45250203));
		L.add(new IntervalSplit(9,45350203,45815521));
		L.add(new IntervalSplit(9,45865521,46216430));
		L.add(new IntervalSplit(9,46266430,46461039));
		L.add(new IntervalSplit(9,46561039,47060133));
		L.add(new IntervalSplit(9,47160133,47317679));
		L.add(new IntervalSplit(9,65467679,65918360));
		L.add(new IntervalSplit(9,65968360,66192215));
		L.add(new IntervalSplit(9,66242215,66404656));
		L.add(new IntervalSplit(9,66454656,66614195));
		L.add(new IntervalSplit(9,66664195,66863343));
		L.add(new IntervalSplit(9,66913343,67107834));
		L.add(new IntervalSplit(9,67207834,67366296));
		L.add(new IntervalSplit(9,67516296,67987998));
		L.add(new IntervalSplit(9,68137998,68514181));
		L.add(new IntervalSplit(9,68664181,68838946));
		L.add(new IntervalSplit(9,68988946,69278385));
		L.add(new IntervalSplit(9,69328385,70010542));
		L.add(new IntervalSplit(9,70060542,70218729));
		L.add(new IntervalSplit(9,70318729,70506535));
		L.add(new IntervalSplit(9,70556535,70735468));
		L.add(new IntervalSplit(9,70835468,92343416));
		L.add(new IntervalSplit(9,92443416,92528796));
		L.add(new IntervalSplit(9,92678796,133073060));
		L.add(new IntervalSplit(9,133223060,137041193));
		L.add(new IntervalSplit(9,137091193,139166997));
		L.add(new IntervalSplit(9,139216997,141153431));
		L.add(new IntervalSplit(10,60001,17974674));
		L.add(new IntervalSplit(10,18024675,38818835));
		L.add(new IntervalSplit(10,38868835,39154935));
		L.add(new IntervalSplit(10,42354935,42546687));
		L.add(new IntervalSplit(10,42596687,46426964));
		L.add(new IntervalSplit(10,46476964,47429169));
		L.add(new IntervalSplit(10,47529169,47792476));
		L.add(new IntervalSplit(10,47892476,48055707));
		L.add(new IntervalSplit(10,48105707,49095536));
		L.add(new IntervalSplit(10,49195536,51137410));
		L.add(new IntervalSplit(10,51187410,51398845));
		L.add(new IntervalSplit(10,51448845,125869472));
		L.add(new IntervalSplit(10,125919472,128616069));
		L.add(new IntervalSplit(10,128766069,133381404));
		L.add(new IntervalSplit(10,133431404,133677527));
		L.add(new IntervalSplit(10,133727527,135524747));
		L.add(new IntervalSplit(11,60001,1162758));
		L.add(new IntervalSplit(11,1212759,50783853));
		L.add(new IntervalSplit(11,51090853,51594205));
		L.add(new IntervalSplit(11,54694205,69089801));
		L.add(new IntervalSplit(11,69139801,69724695));
		L.add(new IntervalSplit(11,69774695,87688378));
		L.add(new IntervalSplit(11,87738378,96287584));
		L.add(new IntervalSplit(11,96437584,134946516));
		L.add(new IntervalSplit(12,60001,95738));
		L.add(new IntervalSplit(12,145739,7189876));
		L.add(new IntervalSplit(12,7239876,34856694));
		L.add(new IntervalSplit(12,37856694,109373470));
		L.add(new IntervalSplit(12,109423470,122530623));
		L.add(new IntervalSplit(12,122580623,132706992));
		L.add(new IntervalSplit(12,132806992,133841895));
		L.add(new IntervalSplit(13,19020001,86760323));
		L.add(new IntervalSplit(13,86910324,112353994));
		L.add(new IntervalSplit(13,112503994,114325993));
		L.add(new IntervalSplit(13,114425993,114639948));
		L.add(new IntervalSplit(13,114739948,115109878));
		L.add(new IntervalSplit(14,19000001,107289539));
		L.add(new IntervalSplit(15,20000001,20894632));
		L.add(new IntervalSplit(15,20935075,21398819));
		L.add(new IntervalSplit(15,21885000,22212114));
		L.add(new IntervalSplit(15,22262114,22596193));
		L.add(new IntervalSplit(15,22646193,23514853));
		L.add(new IntervalSplit(15,23564853,29159443));
		L.add(new IntervalSplit(15,29209443,82829645));
		L.add(new IntervalSplit(15,82879645,84984473));
		L.add(new IntervalSplit(15,85034473,102521392));
		L.add(new IntervalSplit(16,60001,8636920));
		L.add(new IntervalSplit(16,8686921,34023150));
		L.add(new IntervalSplit(16,34173150,35285801));
		L.add(new IntervalSplit(16,46385801,88389383));
		L.add(new IntervalSplit(16,88439383,90294753));
		L.add(new IntervalSplit(17,0,296626));
		L.add(new IntervalSplit(17,396626,21566608));
		L.add(new IntervalSplit(17,21666608,22263006));
		L.add(new IntervalSplit(17,25263006,34675848));
		L.add(new IntervalSplit(17,34725848,62410760));
		L.add(new IntervalSplit(17,62460760,77546461));
		L.add(new IntervalSplit(17,77596461,79709049));
		L.add(new IntervalSplit(17,79759049,81195210));
		L.add(new IntervalSplit(18,10001,15410897));
		L.add(new IntervalSplit(18,18510898,52059136));
		L.add(new IntervalSplit(18,52209136,72283353));
		L.add(new IntervalSplit(18,72333353,75721820));
		L.add(new IntervalSplit(18,75771820,78017248));
		L.add(new IntervalSplit(19,60001,7346003));
		L.add(new IntervalSplit(19,7396004,8687198));
		L.add(new IntervalSplit(19,8737198,20523415));
		L.add(new IntervalSplit(19,20573415,24631782));
		L.add(new IntervalSplit(19,27731782,59118983));
		L.add(new IntervalSplit(20,60001,26319568));
		L.add(new IntervalSplit(20,29419569,29653908));
		L.add(new IntervalSplit(20,29803908,34897085));
		L.add(new IntervalSplit(20,34947085,61091437));
		L.add(new IntervalSplit(20,61141437,61213369));
		L.add(new IntervalSplit(20,61263369,62965520));
		L.add(new IntervalSplit(21,9411194,9595547));
		L.add(new IntervalSplit(21,9645548,9775437));
		L.add(new IntervalSplit(21,9825437,10034920));
		L.add(new IntervalSplit(21,10084920,10215976));
		L.add(new IntervalSplit(21,10365976,10647896));
		L.add(new IntervalSplit(21,10697896,11188129));
		L.add(new IntervalSplit(21,14338129,42955559));
		L.add(new IntervalSplit(21,43005559,44632664));
		L.add(new IntervalSplit(21,44682664,48119895));
		L.add(new IntervalSplit(22,16050001,16697849));
		L.add(new IntervalSplit(22,16847850,20509431));
		L.add(new IntervalSplit(22,20609431,50364777));
		L.add(new IntervalSplit(22,50414777,51244566));
		
		L.add(new IntervalSplit("X",60001,94820));
		L.add(new IntervalSplit("X",144821,231384));
		L.add(new IntervalSplit("X",281384,1047557));
		L.add(new IntervalSplit("X",1097557,1134113));
		L.add(new IntervalSplit("X",1184113,1264234));
		L.add(new IntervalSplit("X",1314234,2068238));
		L.add(new IntervalSplit("X",2118238,7623882));
		L.add(new IntervalSplit("X",7673882,10738674));
		L.add(new IntervalSplit("X",10788674,37098256));
		L.add(new IntervalSplit("X",37148256,49242997));
		L.add(new IntervalSplit("X",49292997,49974173));
		L.add(new IntervalSplit("X",50024173,52395914));
		L.add(new IntervalSplit("X",52445914,58582012));
		L.add(new IntervalSplit("X",61682012,76653692));
		L.add(new IntervalSplit("X",76703692,113517668));
		L.add(new IntervalSplit("X",113567668,115682290));
		L.add(new IntervalSplit("X",115732290,120013235));
		L.add(new IntervalSplit("X",120063235,143507324));
		L.add(new IntervalSplit("X",143557324,148906424));
		L.add(new IntervalSplit("X",148956424,149032062));
		L.add(new IntervalSplit("X",149082062,152277099));
		L.add(new IntervalSplit("X",152327099,155260560));
		L.add(new IntervalSplit("Y",10001,44820));
		L.add(new IntervalSplit("Y",94821,181384));
		L.add(new IntervalSplit("Y",231384,997557));
		L.add(new IntervalSplit("Y",1047557,1084113));
		L.add(new IntervalSplit("Y",1134113,1214234));
		L.add(new IntervalSplit("Y",1264234,2018238));
		L.add(new IntervalSplit("Y",2068238,8914955));
		L.add(new IntervalSplit("Y",8964955,9241322));
		L.add(new IntervalSplit("Y",9291322,10104553));
		L.add(new IntervalSplit("Y",13104553,13143954));
		L.add(new IntervalSplit("Y",13193954,13748578));
		L.add(new IntervalSplit("Y",13798578,20143885));
		L.add(new IntervalSplit("Y",20193885,22369679));
		L.add(new IntervalSplit("Y",22419679,23901428));
		L.add(new IntervalSplit("Y",23951428,28819361));
		L.add(new IntervalSplit("Y",58819361,58917656));
		L.add(new IntervalSplit("Y",58967656,59363566));
		return L;
		}
	
	public static void main(String[] args) {
		new NgsWorkflow().instanceMainWithExit(args);
	}
	
}
