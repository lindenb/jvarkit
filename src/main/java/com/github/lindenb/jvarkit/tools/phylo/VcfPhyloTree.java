/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.phylo;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.ContigPosRef;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.bdb.AlleleBinding;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.CursorConfig;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

public class VcfPhyloTree extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfPhyloTree.class).make();

	private enum GType
		{
		HOM_REF,
		HOM_ALT,
		HET
		//,UNAVAILABLE
		}
	
	
	
	private static class Call
		{
		int sample_index;
		GType gtype;
		int depth;
		int qual;
		}
	
	private static class CallBinding
		extends TupleBinding<Call>
		{
		private final GType[] GTYPES=GType.values();
		@Override
		public Call entryToObject(TupleInput in) {
			int sample_index= in.readInt();
			if(sample_index==-1) return null;
			Call call= new Call();
			call.sample_index=sample_index;
			call.gtype=GTYPES[(int)in.readByte()];
			call.depth = in.readInt();
			call.qual = in.readInt();
			return call;
			}
		@Override
		public void objectToEntry(Call g, TupleOutput w)
			{
			if(g==null || g.sample_index==-1)
				{	
				w.writeInt(-1);
				return;
				}
			w.writeInt(g.sample_index);
			w.writeByte((byte)g.gtype.ordinal());
			w.writeInt(g.depth);
			w.writeInt(g.qual);
			}
		}

	
	private class Calls
		{
		private List<Call> calls=new ArrayList<>(VcfPhyloTree.this.sampleList.size());
		
		
		Call get(int sampleIndex)
			{
			return sampleIndex<calls.size()?calls.get(sampleIndex):null;
			}
		GType getGenotype(int sampleIndex)
			{
			Call c=get(sampleIndex);
			return c==null?GType.HOM_REF:c.gtype;
			}
		void visit(VariantContext ctx)
			{
			for(int i=0;i<  VcfPhyloTree.this.sampleList.size();++i)
				{
				String sampleName= VcfPhyloTree.this.sampleList.get(i);
				Genotype g = ctx.getGenotype(sampleName);
				Call newcall=new Call();
				newcall.sample_index=i;
				if(g.hasDP()) newcall.depth=g.getDP();
				if(g.hasGQ()) newcall.qual=g.getGQ();
				
				switch(g.getType())
					{
					case HET: newcall.gtype =GType.HET;break;
					case HOM_REF: newcall.gtype =GType.HOM_REF;break;
					case HOM_VAR: newcall.gtype =GType.HOM_ALT;break;
					default: newcall=null; break;
					}
				if(newcall==null) continue;
				
				Call old= this.get(i);
				if( old!=null && 
					old.qual!=-1 && newcall.qual!=-1 &&
					old.qual > newcall.qual
					)
					{
					continue;
					}
				while(this.calls.size()<=i) this.calls.add(null);
				this.calls.set(i, newcall);
				}
			}
		}

	private class CallsBinding
	extends TupleBinding<Calls>
		{
		private CallBinding callBinding=new CallBinding();
		@Override
		public Calls entryToObject(TupleInput in)
			{
			Calls calls=new Calls();
			int n=in.readInt();
			for(int i=0;i< n;++i)
				{
				calls.calls.add(this.callBinding.entryToObject(in));
				}
			return calls;
			}
		@Override
		public void objectToEntry(Calls o, TupleOutput w)
			{
			w.writeInt(o.calls.size());
			for(int i=0;i< o.calls.size();++i)
				{
				this.callBinding.objectToEntry(o.calls.get(i),w);
				}
			}
		}
	private  final CallsBinding callsBindingInstance=new CallsBinding();

	

	
	private static class ChromPosRefBinding
		extends TupleBinding<ContigPosRef>
		{
		AlleleBinding alleleBinding=new AlleleBinding();
		@Override
		public ContigPosRef entryToObject(TupleInput in)
			{
			final String c = in.readString();
			int p = in.readInt();
			final Allele a = alleleBinding.entryToObject(in);
			return new ContigPosRef(c,p,a);
			}
		@Override
		public void objectToEntry(ContigPosRef o, TupleOutput out) {
			out.writeString(o.getContig());
			out.writeInt(o.getStart());
			alleleBinding.objectToEntry(o.getReference(), out);
			}
		}
	private ChromPosRefBinding chromPosRefBindingInstance=new ChromPosRefBinding();
	
	
	
	private long ID_GENERATOR=0L;
	
	
	
	private abstract class TreeNode
		{
		long nodeid=(++ID_GENERATOR);
		double weight = 0L;
		private Set<Integer> _cacheSamples=null;
		
		abstract void collectSampleIds(Set<Integer> set);
		public Set<Integer> getSamples()
			{
			if(_cacheSamples==null)
				{
				_cacheSamples=new HashSet<Integer>();
				collectSampleIds(_cacheSamples);
				}
			return _cacheSamples;
			}
		
		public Set<String> getSampleNames()
			{
			Set<String> h=new TreeSet<String>();
			for(Integer idx:this.getSamples()) h.add(VcfPhyloTree.this.sampleList.get(idx));
			return h;
			}
		
		protected abstract Counter<GType> getGenotypeTypes( Counter<GType> counter,final Calls ctx);
		
		private double distance(GType t1,GType t2)
			{	
			final double dAA_AA=0.0;
			final double dAA_AB=5.0;
			final double dAA_BB=15.0;

			
			switch(t1)
				{
				case HET:
					{
					switch(t2)
						{
						case HET: return dAA_AA;
						case HOM_ALT: return dAA_AB;
						case HOM_REF: return dAA_AB;
						default:break;
						}
					break;
					}
				case HOM_ALT:
					{
					switch(t2)
						{
						case HET: return dAA_AB;
						case HOM_ALT: return dAA_AA;
						case HOM_REF: return dAA_BB;
						default:break;
						}
					break;
					}
				case HOM_REF:
					{
					switch(t2)
						{
						case HET: return dAA_AB;
						case HOM_ALT: return dAA_BB;
						case HOM_REF: return dAA_AA;
						default:break;
						}
					break;
					}
				default:break;
				}
			throw new RuntimeException(t1.toString()+" "+t2);
			}
		
		double distance(TreeNode other,final Calls ctx)
			{
			Counter<GType> x1=new Counter<GType>();
			Counter<GType> x2=new Counter<GType>();
			x1= this.getGenotypeTypes(x1,ctx);
			x2= other.getGenotypeTypes(x2,ctx);
			if(x1.isEmpty() && x2.isEmpty()) return 0.0;
			if(x1.isEmpty() || x2.isEmpty()) return 1.0;
			double n=0;
			long m=0;
			for(GType t1:x1.keySet())
				{
				for(GType t2:x2.keySet())
					{
					n+=distance(t1,t2)*x1.count(t1)*x2.count(t2);
					m++;
					}
				}
			
			return n/m;
			/*
			if(x1.keySet().equals(x2.keySet())) return 0L;
			
			GType t1= x1.getMostFrequent();
			GType t2= x2.getMostFrequent();
			return distance(t1,t2);
			*/
			}

		
		boolean containsAnySampleIds(Set<Integer> set)
			{
			for(Integer s: getSamples())
				{
				if(set.contains(s)) return true;
				}
			return false;
			}
		abstract void writeGraphizDot(PrintStream out);
		
		abstract void writeGexfNodes(XMLStreamWriter out) throws XMLStreamException;
		abstract void writeGexfEdges(XMLStreamWriter out) throws XMLStreamException;
		
		abstract void writeNewick(PrintStream out);

		
		public String getNodeId()
			{
			return "node"+this.nodeid;
			}
		abstract public String toString();
		
		}
	
	private class OneSample extends TreeNode
		{
		private int sample_index;
		OneSample(int sample_index)
			{
			this.sample_index = sample_index;
			}
		String getSampleName()
			{
			return VcfPhyloTree.this.sampleList.get(this.sample_index);
			}
		
		@Override
		protected Counter<GType> getGenotypeTypes(Counter<GType> c,final Calls calls)
			{
			c.incr(calls.getGenotype(this.sample_index));
			return c;
			}
		
		@Override
		void collectSampleIds(Set<Integer> set)
			{
			set.add(this.sample_index);
			}
		
		@Override
		void writeGraphizDot(PrintStream out) {
			out.println(getNodeId()+"[label=\""+this.getSampleName()+"\"];");
			}
		@Override
		 void writeGexfNodes(XMLStreamWriter w) throws XMLStreamException
			{
			w.writeEmptyElement("node");
				w.writeAttribute("id",getNodeId());
				w.writeAttribute("label",this.getSampleName());
			}
		
		@Override
		void writeGexfEdges(XMLStreamWriter w) throws XMLStreamException
			{
			//empty
			}
		@Override
		void writeNewick(PrintStream out)
			{
			out.print(this.getSampleName());
			}
		@Override
		public String toString() {
			return this.getSampleName()+" ("+this.weight+")";
			}
		}
	
	

	
	private class MergedNodes extends TreeNode
		{
		private TreeNode t1;
		private TreeNode t2;
		MergedNodes(TreeNode t1,TreeNode t2)
			{
			this.t1=t1;
			this.t2=t2;
			}
		
		@Override
		void collectSampleIds(Set<Integer> set) {
			t1.collectSampleIds(set);
			t2.collectSampleIds(set);
			}
		
		
		@Override
		protected Counter<GType> getGenotypeTypes(Counter<GType> c,final Calls ctx)
			{
			t1.getGenotypeTypes(c,ctx);
			t2.getGenotypeTypes(c,ctx);
			
			return c;
			}

		@Override
		void writeGraphizDot(PrintStream out)
			{
			t1.writeGraphizDot(out);
			t2.writeGraphizDot(out);
			out.println(this.getNodeId()+"[shape=point];");
			String label="[weight="+(1+(int)this.weight)+",label=\""+(int)this.weight+"\"]";
			out.println(t1.getNodeId()+" -- "+ this.getNodeId()+label+";");
			out.println(t2.getNodeId()+" -- "+ this.getNodeId()+label+";");
			
			}
		@Override
		 void writeGexfNodes(XMLStreamWriter w) throws XMLStreamException
			{
			w.writeEmptyElement("node");
				w.writeAttribute("id",getNodeId());
				w.writeAttribute("label","");
			this.t1.writeGexfNodes(w);
			this.t2.writeGexfNodes(w);
			}
		
		@Override
		void writeGexfEdges(XMLStreamWriter w) throws XMLStreamException
			{
			for(int i=0;i<2;++i)
				{
				TreeNode c=(i==0?t1:t2);
				String nid= c.getNodeId();
				w.writeStartElement("edge");
				w.writeAttribute("id","E"+this.nodeid+"_"+nid);
				w.writeAttribute("type","directed");
				w.writeAttribute("source",nid);
				w.writeAttribute("target",this.getNodeId());
				w.writeAttribute("weight",String.valueOf((int)this.weight));
				/*
				w.writeStartElement("attvalues");
				w.writeEmptyElement("attvalue");
				 w.writeAttribute("for","weight");
				w.writeAttribute("value",String.valueOf(this.score));
				w.writeEndElement();//attvalues
				*/
				w.writeEndElement();
				}
			this.t1.writeGexfEdges(w);
			this.t2.writeGexfEdges(w);
			}
		@Override
		void writeNewick(PrintStream w)
			{
			w.print("(");
			this.t1.writeNewick(w);
			w.print(",");
			this.t2.writeNewick(w);
			w.print(")");
			}
		
		@Override
		public String toString()
			{
			Set<String> s=getSampleNames();
			if(s.size()>10) return ""+s.size()+" samples ("+this.weight+")";
			return Arrays.toString(s.toArray())+"("+this.weight+")";
			}
		}
	

	
	/** map TidPos to CovRow */
	private Environment env;
	/** map TidPos to CovRow */
	private Database pos2cov;
	/** all samples */
	private Map<String,Integer> sample2col=new HashMap<String,Integer>();
	private List<String> sampleList =new ArrayList<String>();
	
	/** count snps */
	private long snp_count=0L;
	
	private VcfPhyloTree()
		{
		
		}
	
	private void cleanup(Transaction txn)
		{
		long count_deleted=0;
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		Cursor cursor=null;
		cursor=this.pos2cov.openCursor(txn, null);
		while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
			{
			/*
			Calls row = this.callsBindingInstance.entryToObject(data);
			
			boolean mixed=false;
			for(String sample:this.sampleList)
				{
				Genotype g=row.get(sample);
				
				if(g!=null && g.isMixed())
					{
					mixed=true;
					break;
					}
				}
			//remove mixed calls
			if(mixed)
				{
				cursor.delete();
				++count_deleted;
				continue;
				}*/
			}
		cursor.close();
		LOG.info("After cleanup, removed "+count_deleted+" snps.");
		}
	private final int MAX_NODES=2000000;
	
	private TreeNode matrix(Transaction txn, List<TreeNode> parents)
		{
		long last=System.currentTimeMillis();
		
		/* faster computing ? , parse BerkeleyDB one time */
		if(Math.pow(parents.size(),2)/2.0 < MAX_NODES)
			{
			CursorConfig cursorConfig=CursorConfig.DEFAULT;
			List<MergedNodes> nodes=new ArrayList<>();
			for(int x=0;x< parents.size();++x)
				{
				TreeNode tn1= parents.get(x);
				for(int y=x+1;y< parents.size();++y)
					{
					TreeNode tn2= parents.get(y);
					nodes.add(new MergedNodes(tn1, tn2));
					}
				}
			DatabaseEntry key=new DatabaseEntry();
			DatabaseEntry data=new DatabaseEntry();
			Cursor cursor=this.pos2cov.openCursor(txn,cursorConfig);
			long nVar=0L;
			while(cursor.getNext(key, data,LockMode.DEFAULT )==OperationStatus.SUCCESS)
				{
				long now=System.currentTimeMillis();
				if(now-last > 10*1000)
					{
					LOG.info("iteration "+nVar+"/"+this.snp_count+" for "+parents.size());
					last=now;
					}
				
				Calls row = this.callsBindingInstance.entryToObject(data);
				for(MergedNodes merged : nodes)
					{
					merged.weight += merged.t1.distance(merged.t2,row);
					}
				++nVar;
				}
			cursor.close();
			
			Collections.sort(nodes,new Comparator<MergedNodes>()
				{
				public int compare(MergedNodes o1, MergedNodes o2)
					{
					if(o1.weight< o2.weight) return -1;
					if(o1.weight> o2.weight) return 1;
					return 0;
					}
				});
			
			if(nodes.size()>1 &&
					nodes.get(0).weight == nodes.get(1).weight
					)
				{
				LOG.warning("Score ambiguity for "+nodes.get(0)+" and "+nodes.get(1));
				}
				
			return nodes.get(0);
			}
		else
			{
			Cursor cursor=this.pos2cov.openCursor(txn, null);;
			MergedNodes best=null;
			for(int x=0;x< parents.size();++x)
				{
				TreeNode tn1= parents.get(x);
				
				
				for(int y=x+1;y< parents.size();++y)
					{
					long now=System.currentTimeMillis();
					if(now-last > 10*1000)
						{
						LOG.info("iteration ["+(x+1)+","+(y+1)+"]/"+parents.size());
						last=now;
						}
					
					TreeNode tn2= parents.get(y);
					DatabaseEntry key=new DatabaseEntry();
					DatabaseEntry data=new DatabaseEntry();
					MergedNodes merged= new MergedNodes(tn1, tn2);
					boolean first=true;
					
					
					while((first?cursor.getFirst(key, data, LockMode.DEFAULT):cursor.getNext(key, data, LockMode.DEFAULT))==OperationStatus.SUCCESS)
						{
						first=false;
						Calls row = this.callsBindingInstance.entryToObject(data);
						
						merged.weight += tn1.distance(tn2,row);
						}
					
					if(best==null || best.weight > merged.weight)
						{
						best=merged;
						}
					}
				}
			cursor.close();
			return best;
			}
		
		}
	
	private TreeNode matrix(Transaction txn)
		{
		ArrayList<TreeNode> nodes=new ArrayList<TreeNode>();
		//initialize population
		for(int sample_index=0; sample_index< this.sampleList.size();++sample_index)
			{
			nodes.add(new OneSample(sample_index));
			}
		
		while(nodes.size()>1)
			{
			long now= System.currentTimeMillis();
			LOG.info("nodes.count= "+nodes.size());
			
			/* no, want to get the distance between the two nodes 
			if(nodes.size()==2)
				{
				return new MergedNodes(
						nodes.get(0),
						nodes.get(1)
						);
				}
			*/
			
			
			TreeNode best = matrix(txn, nodes);
			LOG.info(best);
			Set<Integer> bestsamples= best.getSamples();
			
			ArrayList<TreeNode> newnodes=new ArrayList<TreeNode>();

			for(TreeNode tn:nodes)
				{
				if(tn.containsAnySampleIds(bestsamples)) continue;
				newnodes.add(tn);
				}
			newnodes.add(best);
			
			nodes=newnodes;
			
			LOG.info("That took :"+ (System.currentTimeMillis()-now)/1000f  +" seconds");
			}
		return nodes.get(0);
		}
	
	
	private void readvcf(Transaction txn,VcfIterator in)
		{
		Cursor cursor=null;
		VCFHeader header = in.getHeader();
		if(!header.hasGenotypingData()) return;
		List<String> sampleNames = header.getSampleNamesInOrder();
		for(String sample: sampleNames)
			{
			if(this.sample2col.containsKey(sample)) continue;
			this.sample2col.put(sample, this.sampleList.size());
			this.sampleList.add(sample);
			}
		
		cursor= this.pos2cov.openCursor(txn, null);

		DatabaseEntry key = new DatabaseEntry();
		DatabaseEntry data = new DatabaseEntry();
		
		while(in.hasNext())
			{
			VariantContext ctx = in.next();
			ContigPosRef chromPosRef=new ContigPosRef(ctx);
			
			
			chromPosRefBindingInstance.objectToEntry(chromPosRef, key);
			
			if(cursor.getSearchKey(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				Calls row = this.callsBindingInstance.entryToObject(data);
				row.visit(ctx);
				this.callsBindingInstance.objectToEntry(row, data);
				cursor.putCurrent(data);
				}
			else
				{
				Calls row=new Calls();
				row.visit(ctx);
				this.callsBindingInstance.objectToEntry(row, data);
				cursor.put(key,data);
				}
			
			}
		cursor.close();
		}
	

	
	
	private enum OUT_FMT
		{
		dot,gexf,newick
		};
	
		
	@Parameter(names="-B",description="bdb home")
	private 	File bdbHome=null;;
	@Parameter(names="-f",description="output format")
	private 	OUT_FMT outformat=OUT_FMT.dot;

	@Override
	public int doWork(List<String> args) {

		if(bdbHome==null)
			{
			LOG.error("undefined bdb home");
			return -1;
			}
		
		
		Transaction txn=null;
		VcfIterator iter=null;
		try
			{
			LOG.info("Opening "+bdbHome);
			EnvironmentConfig cfg=new EnvironmentConfig();
			cfg.setAllowCreate(true);
			cfg.setReadOnly(false);
			//cfg.setTransactional(false);
			
			this.env= new Environment(bdbHome, cfg);
			
			DatabaseConfig cfg2=new DatabaseConfig();
			cfg2.setAllowCreate(true);
			cfg2.setTemporary(true);
			cfg2.setReadOnly(false);
			
			//cfg2.setTransactional(false);
			this.pos2cov=this.env.openDatabase(txn, "vcf", cfg2);
			
			
			if(args.isEmpty())
				{
				LOG.info("reading from stdin");
				iter=super.openVcfIterator(null);
				readvcf(txn,iter);
				iter.close();
				}
			else
				{
				for(String filename:args)
					{
					LOG.info("reading "+filename);
					iter=VCFUtils.createVcfIterator(filename);
					readvcf(txn,iter);
					iter.close();
					}
				}
			
			if(this.sampleList.size()<2)
				{
				LOG.error("Not enoug samples");
				return -1;
				}
			
			
			
			cleanup(txn);
			
			this.snp_count= this.pos2cov.count();
			LOG.info("Samples: "+this.sampleList.size());
			LOG.info("SNPs: "+this.snp_count);

			if(this.snp_count==0L)
				{
				LOG.error("Not enough SNPS.");
				return -1;
				}
			
			TreeNode t = matrix(txn);
			
			switch(outformat)
				{
				case newick:
					{
					t.writeNewick(System.out);
					System.out.print(";");
					break;
					}
				case gexf:
					{
					XMLOutputFactory xof=XMLOutputFactory.newFactory();
					XMLStreamWriter w=xof.createXMLStreamWriter(System.out, "UTF-8");
					w.writeStartDocument("UTF-8", "1.0");
					w.writeStartElement("gexf");
					w.writeDefaultNamespace("http://www.gexf.net/1.2draft");
					w.writeAttribute("version", "1.2");
					w.writeStartElement("meta");
					  w.writeStartElement("creator");
					  w.writeCharacters(getProgramName()+" by Pierre Lindenbaum");
					  w.writeEndElement();
					 
					  w.writeStartElement("description");
					  w.writeCharacters(getProgramCommandLine());
					  w.writeEndElement();
					
					w.writeEndElement();//meta
					  
					  w.writeStartElement("graph");
					  w.writeAttribute("mode", "static");
					  w.writeAttribute("defaultedgetype", "directed");
					  
					  w.writeStartElement("attributes");
					  w.writeAttribute("class", "node");
					  w.writeAttribute("mode", "static");
					  w.writeEndElement();//attributes
						
					  w.writeStartElement("attributes");                                                                                     
					  w.writeAttribute("class", "edge");
					  w.writeAttribute("mode", "static");
    					  /*
                        w.writeEmptyElement("attribute");
						w.writeAttribute("id", "weight");
						w.writeAttribute("title", "Weight");
						w.writeAttribute("type", "float");*/
                      w.writeEndElement();//attributes
					  

					  
					  w.writeStartElement("nodes");
					  t.writeGexfNodes(w);
					  w.writeEndElement();//nodes
					
					  
					  w.writeStartElement("edges");
					  t.writeGexfEdges(w);
					  w.writeEndElement();//edges
					  
					  w.writeEndElement();//graph

					
					w.writeEndElement();//gexf
					w.writeEndDocument();
					w.flush();
					System.out.flush();
					break;
					}
				case dot:
				default:
					{
					System.out.println("graph G {");
					t.writeGraphizDot(System.out);
					System.out.println("}");
					break;
					}
				}
			return 0;
			}
		catch(Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(this.pos2cov);
			CloserUtil.close(this.env);
			}
		}
	public static void main(String[] args) {
		new VcfPhyloTree().instanceMainWithExit(args);
		}
	}
