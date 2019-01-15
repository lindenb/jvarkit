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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC




### Example



```

$ java -Xmx4g -cp ${classpath} \
     org.broadinstitute.gatk.queue.QCommandLine \
     -S ${SV_DIR}/qscript/SVPreprocess.q \
     -S ${SV_DIR}/qscript/SVQScript.q \
     -cp ${classpath} \
     -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
     -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
     -R ${REF} \
     -I bam.list \
     -md output_metadata_directory \
     -bamFilesAreDisjoint true \
     -jobLogDir logDir 2> shell.txt
   
$ java -jar dist/queue2make.jar shell.txt   > shell.mk

```






END_DOC
*/
@Program(name="queue2make",description="Convert Broad/Queue genomestrip Log stream to Makefile.")
public class QueueToMake extends Launcher
	{
	private static final Logger LOG = Logger.build(QueueToMake.class).make();
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	
	private Map<String,Target> file2target = new HashMap<>();
	final List<Target> orderedTargets = new ArrayList<>();

	
	private class Command extends AbstractList<String>
		{
		private List<String> args = new ArrayList<>();
		Command() {}
		Command(final String s)
			{
			int i=0;
			while(i< s.length())
				{
				if(Character.isWhitespace(s.charAt(i))) {
					++i;
					continue;
					}
				else if(s.charAt(i)=='\'') {
					i++;
					final StringBuilder sb =new StringBuilder();
					for(;;) {
						if(i==s.length()) throw new IllegalArgumentException(s);
						if(s.charAt(i)=='\'') break;
						sb.append(s.charAt(i));
						i++;
						}
					i++;//last apos
					this.args.add(sb.toString());
					continue;
					}
				else
					{
					final StringBuilder sb =new StringBuilder();
					while(i< s.length() && !Character.isWhitespace(s.charAt(i)))
						{
						sb.append(s.charAt(i));
						i++;
						}
					this.args.add(sb.toString());
					}
				}
			for(i=0;i< this.args.size();++i) {
				if(args.get(i).startsWith("-Djava.io.tmpdir=")) {
					args.set(i,"-Djava.io.tmpdir=$(dir $@)");
				}
			}
			}
		@Override
		public int size() {
			return args.size();
			}
		@Override
		public String get(int index) {
			return args.get(index);
			}
		
		 void prepend(String...atts) {
	        for(int i=atts.length-1;i>=0;--i) {
	        	this.add(0, atts[i]);
	        	}
      		}
		
		@Override
		public String toString() {
			return String.join(" ", args);
			}
		String last() {
			return this.get(this.size()-1);
		}
		public List<String> allValuesFor(String opt)
	        {
	        final List<String> l= new ArrayList<>();
	        for(int i=0;i+1< this.size();++i)
	            {
	            if(!this.get(i).equals(opt)) continue;
	            i++;
	            String v=this.get(i);
	            if(v.startsWith("-"))   throw new IllegalStateException("Value starting with hyphen in "+this.toString());
	            l.add(v);
	            }
	        return l;
	        }
		@Override
		public void add(int index, String element) {
			this.args.add(index, element);
			}
    
		public  Optional<String> optionalValueFor(final String opt)
	        {
	        final List<String> v= allValuesFor(opt);
	        if(v.isEmpty()) return Optional.empty();
	        if(v.size()==1) return Optional.of(v.get(0));
	        throw new IllegalStateException("found more that one value for "+opt+" in "+this.toString());
	        }
		public String requiredValueFor(final String opt)
	        {
	        final Optional<String> v = optionalValueFor(opt);
	        if(v.isPresent()) return v.get();
	        throw new IllegalStateException("found no value for "+opt+" in "+this.toString());
	        }
		
		@SuppressWarnings("unused")
		private void debug() {
			for(int x=0;x< args.size();++x)
				{
				if(!args.get(x).startsWith("-")) continue;
				if(args.get(x+1).startsWith("-")) continue;
				String prefix=file2target.containsKey(args.get(x+1))?"PREREQUISTE:":"";
				System.err.println(prefix+" "+args.get(x)+" = "+args.get(x+1));
				}
			}
		}
	
    private  class Target
        {
    	private boolean printed=false;
        private final String file;
        private final Command command=new Command();
        private final Set<Target> prerequisites = new HashSet<>();
        Target(String file) {
            this.file=file;
            }
        Target(String file,final Command cmd) {
            this(file);
            this.command.addAll(cmd);
            if(!cmd.isEmpty() && cmd.get(0).equals("java") && cmd.contains("-md") && !cmd.contains("org.broadinstitute.sv.apps.CreateMetaDataDirectory"))
            	{
                this.prerequisites.add( QueueToMake.this.mustExistsTarget(new File(cmd.requiredValueFor("-md"),"success.txt").getPath()));
            	}
            }
        
      
        
        void make(PrintWriter out) {
        	if(printed) return;
        	if(this.prerequisites.isEmpty() && this.command.isEmpty()) {
        		printed=true;
        		return;
        	}
        	String firstPrereq = null;
        	out.print(this.file);
        	out.print(" : ");
        	for(Target dep: this.prerequisites) {
        		if(firstPrereq==null) firstPrereq=dep.file;
        		out.print(" \\\n\t"+dep.file);
        	}
        	if(!this.command.isEmpty()) {
	        	out.println();
	        	out.print("\tmkdir -p $(dir $@) &&");
	        	for(String part: this.command) 
	        		{
	        		out.print(" ");
	        		if(firstPrereq!=null && firstPrereq.equals(part)) {
	        			out.print("$<");
	        		} else if(this.file.equals(part)) {
	        			out.print("$@");
	        		} else {
	        			out.print(part);
	        		}
	        		}
	        	}
        	out.println();
        	printed=true;
        	for(Target t:prerequisites) {
        		t.make(out);
        		}
        	}
        
        Target autofill()
        	{
        	for(int x=0; x +1< this.command.size();++x)
        		{
        		if(!this.command.get(x).startsWith("-")) continue;
        		++x;
        		final String v= this.command.get(x);
        		
        		if(v.startsWith("-")) continue;
        		final Target dep=QueueToMake.this.file2target.get(v);
        		if(dep==this || dep==null) continue;
        		this.prerequisites.add(dep);
        		
        		}
        	for(int x=0; x +1< this.command.size();++x)
	    		{
	    		if(!this.command.get(x).startsWith("-")) continue;
	    		++x;
	    		if(this.command.get(x).startsWith("-")) continue;
	    		String v= this.command.get(x);
	    		if(v.endsWith(".bam")) {
		    		final Target bam=QueueToMake.this.file2target.get(v);
		    		if(bam==null) {
		    			//file not in target, must be a USER initial bam
		    			this.prerequisites.add(QueueToMake.this.mustNotExistsTarget(v, new Command()));
		    			}
		    		v+=".bai";
		    		final Target dep=QueueToMake.this.file2target.get(v);
	        		if(dep==this || dep==null) continue;
	        		this.prerequisites.add(dep);
	    			}
	    		else if(v.endsWith(".gz")) {
		    		v+=".tbi";
		    		final Target dep=QueueToMake.this.file2target.get(v);
	        		if(dep==this || dep==null) continue;
	        		this.prerequisites.add(dep);
	    			}
	    		else if(v.endsWith(".rccache.bin")) {
		    		v+=".idx";
		    		final Target dep=QueueToMake.this.file2target.get(v);
	        		if(dep==this || dep==null) continue;
	        		this.prerequisites.add(dep);
	    			}
	    		}
        	return this;
        	}
        
        @Override
        public int hashCode() { return this.file.hashCode();}
        @Override
        public boolean equals(final Object o) {
            if(o==this) return true;
            if(o==null || !(o instanceof Target)) return false;
            return this.file.equals(Target.class.cast(o).file);
            }
        @Override
        public String toString() {
            return file;
            }
        }
       
    private Target makeTarget(String file)
        {
        Target t= this.file2target.get(file);
        if(t==null) {
            t=new Target(file);
            this.file2target.put(file,t);
            this.orderedTargets.add(t);
            }
        return t;
        }
    
    private Target mustExistsTarget(final String file)
	    {
	    final Target t= this.file2target.get(file);
	    if(t==null) throw new IllegalStateException("Can't find target that should exists:"+file);
	    return t;
	    }
    private Target mustNotExistsTarget(String file,final Command command)
	    {
	    if(this.file2target.containsKey(file)) throw new IllegalStateException("Found target that should NOT exists:"+file);
	    final Target t=new Target(file,command);
	    this.file2target.put(file,t);
        this.orderedTargets.add(t);
        return t;
	    }
   
    private Target mustExistsTargetEndsWith(final String endOfString)
	    {
    	Target t2=null;
    	for(final Target t : this.file2target.values()) {
    		if(t.file.endsWith(endOfString)) 
    			{
    			if(t2!=null)  throw new IllegalStateException("Both "+t+" and "+t2 +" end with "+endOfString);
    			t2=t;
    			}
    		}
	    if(t2==null) throw new IllegalStateException("Can't find target that should exists:"+endOfString);
	    return t2;
	    }
    
    private abstract class CommandDecoder
        {
        abstract boolean canDecode(final Command args);
        abstract void decode(final Command args);
        }
    
    
   
    private final CommandDecoder decoders[];
    public QueueToMake() {
    	this.decoders = new CommandDecoder[]{
    			new CommandDecoder() {
    				@Override
    				boolean canDecode(final Command args) {
    					return args.contains("-md") &&
    							args.contains("org.broadinstitute.sv.apps.CreateMetaDataDirectory");
    					}
    				@Override
    				void decode(final Command args) {
    					File md = new File(args.requiredValueFor("-md"),"success.txt");
    					Target t= mustNotExistsTarget(md.getPath(),args);
    					t.command.add("&&");
    					t.command.add("touch");
    					t.command.add(t.file);
    					}
    				},
    			new CommandDecoder() {
    				@Override
    				boolean canDecode(final Command args) {
    					return args.contains("-I") && 
    						   args.contains("-O") &&
    						   args.contains("-R") &&
    						   args.contains("org.broadinstitute.sv.apps.ExtractBAMSubset");
    					}
    				@Override
    				void decode(final Command args) {
    					Target t = mustNotExistsTarget(args.requiredValueFor("-O"),args);
    					t.prerequisites.addAll(args.allValuesFor("-I").stream().map(T->makeTarget(T)).collect(Collectors.toSet()));
    					t.command.add("&&");
    					t.command.add("samtools");
    					t.command.add("index");
    					t.command.add(t.file);
    					}
    				},
    			new CommandDecoder() {
    				@Override
    				boolean canDecode(final Command args) {
    					return args.contains("-O") &&
    						   args.contains("-R") &&
    						   args.contains("org.broadinstitute.sv.apps.ComputeGenomeSizes");
    					}
    				@Override
    				void decode(final Command args) {
    					mustNotExistsTarget(args.requiredValueFor("-O"),args);
    					}
    				},
    			new CommandDecoder() {
    					@Override
    					boolean canDecode(final Command args) {
    						return args.contains("org.broadinstitute.sv.main.SVCommandLine") &&
    							   args.contains("-T") &&
    							   args.contains("-md") &&
    							   args.contains("-depthFile") &&
    							   args.contains("-spanFile") &&
    							   args.contains("-gcProfileFile") &&
    							   args.contains("-readCountFile") &&
    							   args.requiredValueFor("-T").equals("ComputeMetadataWalker")
    								;
    						}
    					@Override
    					void decode(final Command args) {
    						for(Target t:file2target.values()) {
    							LOG.info("vailable:\t"+t);
    							}
    						Target depth = mustNotExistsTarget(args.requiredValueFor("-depthFile"), args).autofill();
    						depth.command.prepend("mkdir","-p",
										"$(dir "+args.requiredValueFor("-depthFile")+")",
										"$(dir "+args.requiredValueFor("-spanFile")+")",
    									"$(dir "+args.requiredValueFor("-gcProfileFile")+")",
    									"$(dir "+args.requiredValueFor("-readCountFile")+")",
    									"&&"
    									);
    						
    						File bamFile= new File(args.requiredValueFor("-I"));
    						
    						
    						depth.prerequisites.add(mustExistsTarget(
    								args.requiredValueFor("-md")+"/isd.stats.dat"
    								));
    						
    						depth.prerequisites.add(mustExistsTarget(
    								args.requiredValueFor("-md")+"/isd/"+
    								bamFile.getName().replaceAll(".bam$", ".dist.bin")
    								));
    						
    						Target span = mustNotExistsTarget(args.requiredValueFor("-spanFile"),
    								new Command("touch -c "+args.requiredValueFor("-spanFile")));
    						span.prerequisites.add(depth);
    						Target gc = mustNotExistsTarget(args.requiredValueFor("-gcProfileFile"),
    								new Command("touch -c "+args.requiredValueFor("-gcProfileFile")));
    						gc.prerequisites.add(depth);
    						Target rc = mustNotExistsTarget(args.requiredValueFor("-readCountFile"),
    								new Command("touch -c "+args.requiredValueFor("-readCountFile")));
    						rc.prerequisites.add(depth);
    						
    						
    						}
    					},
    			new CommandDecoder() {
    						@Override
    						boolean canDecode(final Command args) {
    							return args.contains("-O") &&
    								   args.contains("org.broadinstitute.sv.apps.ComputeSampleReadDepth")
    								  ;
    							}
    						@Override
    						void decode(final Command args) {
    							Target t= mustNotExistsTarget(args.requiredValueFor("-O"),args);
    							t.autofill();
    							String profileDirectory=args.requiredValueFor("-profileDirectory");
    						
    							
    							for(String seq:args.allValuesFor("-sequence")) {
    								for(Target t2: file2target.values()) {
    									if(!t2.file.contains(profileDirectory)) continue;
    									if(!t2.file.contains("/profile_seq_"+seq+"_")) continue;
    									if(!t2.file.endsWith(".dat.gz")) continue;
    									Target t3= mustExistsTarget(t2.file+".tbi");
    									t.prerequisites.add(t3);
    									}
    								
    								
    								}
    							}
    						},
    			new CommandDecoder() {
    							@Override
    							boolean canDecode(final Command args) {
    								return args.contains("-O") &&
    										args.contains("org.broadinstitute.sv.apps.ComputeDepthProfiles")
    									   ;
    								}
    							@Override
    							void decode(final Command args) {
    								final Target t=mustNotExistsTarget(args.requiredValueFor("-O"),args).autofill();
    								Target rccacheIdx=null;
    								for(final Target t2:QueueToMake.this.file2target.values()) {
    									if(t2.file.endsWith("/rccache.bin.idx")) {
    										if(rccacheIdx!=null){
    											throw new IllegalArgumentException("multiple "+t2+" and "+rccacheIdx);
    										}
    										rccacheIdx=t2;
    									}
    									}
									if(rccacheIdx==null){
										throw new IllegalArgumentException("cannot find rccache.bin.idx" );
									}
									t.prerequisites.add(rccacheIdx);
    								}
    							},
    			new CommandDecoder() {
					@Override
					boolean canDecode(final Command args) {
						return args.contains("-O") &&
							   args.contains("org.broadinstitute.sv.apps.CallSampleGender")
							   ;
						}
					@Override
					void decode(final Command args) {
						final Target t=mustNotExistsTarget(args.requiredValueFor("-O"),args).autofill();
						
						t.prerequisites.add(QueueToMake.this.mustExistsTargetEndsWith("/rccache.bin.idx"));
						t.prerequisites.add(QueueToMake.this.mustExistsTargetEndsWith("/gcprofiles.zip"));
						t.prerequisites.add(QueueToMake.this.mustExistsTargetEndsWith("/depth.dat"));
						t.prerequisites.add(QueueToMake.this.mustExistsTargetEndsWith("/spans.dat"));
						
						}
					},
    			new CommandDecoder() {
					@Override
					boolean canDecode(final Command args) {
						return args.contains("-O") &&
							   (
							   args.contains("org.broadinstitute.sv.apps.ComputeGCProfiles") || 
							   args.contains("org.broadinstitute.sv.apps.ReduceInsertSizeHistograms") ||
							   args.contains("org.broadinstitute.sv.apps.MergeInsertSizeDistributions") ||
							   args.contains("org.broadinstitute.sv.apps.ComputeInsertStatistics") ||
							   args.contains("org.broadinstitute.sv.apps.IndexReadCountFile") ||
							   args.contains("org.broadinstitute.sv.apps.MergeReadDepthCoverage") ||
							   args.contains("org.broadinstitute.sv.apps.MergeReadSpanCoverage") ||
							   args.contains("org.broadinstitute.sv.apps.MergeGCProfiles") ||
							   args.contains("org.broadinstitute.sv.queue.WriteFileList") ||
							   args.contains("org.broadinstitute.sv.apps.MergeReadCounts") || 
							   args.contains("org.broadinstitute.sv.apps.ComputeDepthProfiles") || 
							   args.contains("org.broadinstitute.sv.apps.ComputeDepthProfiles") || 
							   ( args.contains("org.broadinstitute.sv.main.SVCommandLine") && args.contains("-T") &&
									  (
									  args.requiredValueFor("-T").equals("ComputeInsertSizeHistogramsWalker") 
									  ))
							   );
						}
					@Override
					void decode(final Command args) {
						mustNotExistsTarget(args.requiredValueFor("-O"),args).autofill();
						}
					},
    			new CommandDecoder() {
            				@Override
            				boolean canDecode(final Command args) {
            					return args.size()>2 && 
            						   args.get(0).equals("samtools") &&
            						   args.get(1).equals("index")
            						   ;
            					}
            				@Override
            				void decode(final Command args) {
            					Target t = mustNotExistsTarget(args.get(2)+".bai",args);
            					t.prerequisites.add(mustExistsTarget(args.get(2)));
            					}
            				},
    			new CommandDecoder() {
    				@Override
    				boolean canDecode(final Command args) {
    					return args.size()>2 && 
    						   args.get(0).equals("tabix") &&
    						   args.last().endsWith(".gz")
    						   ;
    					}
    				@Override
    				void decode(final Command args) {
    					Target t = mustNotExistsTarget(args.last()+".tbi",args);
    					t.prerequisites.add(mustExistsTarget(args.last()));
    					}
    				},
    			new CommandDecoder() {
        				@Override
        				boolean canDecode(final Command args) {
        					return args.size()>3 && 
        						   args.get(0).equals("Rscript") &&
        						   args.get(1).endsWith("plot_chr_vs_chr_readdepth.R")
        						   ;
        					}
        				@Override
        				void decode(final Command args) {
        					Target t = mustNotExistsTarget(args.get(3),args);
        					t.prerequisites.add(mustExistsTarget(args.get(2)));
        					}
        				}
    	        };
    	}
    
    @Override
    public int doWork(List<String> args) {
    BufferedReader in=null;
    	PrintWriter out=null;
    	try {
    		final String inputName = oneFileOrNull(args);
    		final String pendingLine=" QGraph - Pending:";
			in = (inputName==null?
					IOUtils.openStreamForBufferedReader(stdin()):
					IOUtils.openURIForBufferedReading(inputName)
					);
			String line;
			while(((line=in.readLine()))!=null)
				{
				if(!line.startsWith("INFO")) continue;
				int pending= line.indexOf(pendingLine);
				if(pending==-1) continue;
				line= line.substring(pending+pendingLine.length()).trim();
				final Command command = new Command(line);
				int x=0;
				for(x=0;x<this.decoders.length;++x)
					{
					if(this.decoders[x].canDecode(command)) {
						this.decoders[x].decode(command);
						break;
						}
					}
				if(x==this.decoders.length) {
					LOG.error("Cannot decode:\n"+command);
					return -1;
					}
				}
			in.close();in=null;
			
			out = super.openFileOrStdoutAsPrintWriter(outputFile);
			
			out.println(".PHONY:all");
			out.print("all:");
			for(final Target t1:this.file2target.values()) {
				boolean top=true;
				for(final Target t2:this.file2target.values()) {
						if(t2.prerequisites.contains(t1)) {
						top=false;
						break;
						}
					}
				if(top) {
					out.print(" \\\n\t"+t1.file);
				}
				}
			out.println();
			
			for(final Target t:this.orderedTargets) {
				t.make(out);
			}
			out.flush();
			out.close();
			LOG.info("done");
			return RETURN_OK;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
			}
    	finally
    		{
    		CloserUtil.close(in);
    		}
    	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new QueueToMake().instanceMainWithExit(args);

	}

}
