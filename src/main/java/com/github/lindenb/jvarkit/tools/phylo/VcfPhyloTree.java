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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

public class VcfPhyloTree extends AbstractCommandLineProgram
	{
	/** position on genome */
	private static class ChromPosRef
		{
		String chrom;
		int pos1;
		Allele ref;
		}
	
	
	private static class AlleleBinding extends TupleBinding<Allele>
		{
		@Override
		public Allele entryToObject(TupleInput in) {
			int n= in.readInt();
			byte bases[]=new byte[n];
			in.read(bases);
			boolean isRef=in.readBoolean();
			return Allele.create(bases, isRef);
			}
		@Override
		public void objectToEntry(Allele a, TupleOutput w)
			{
			byte bases[]=a.getDisplayBases();
			if(bases.length==0) throw new RuntimeException("n=0 for "+a);
			w.writeInt(bases.length);
			w.write(bases);
			w.writeBoolean(a.isReference());
			}
		}

	
	private static class GenotypeBinding extends TupleBinding<Genotype>
		{
		private AlleleBinding alleleBinding=new AlleleBinding();
		
		@Override
		public Genotype entryToObject(TupleInput in) {
			String sampleName = in.readString();
			int n=in.readInt();
			List<Allele> alleles=new ArrayList<Allele>(n);
			for(int i=0;i< n;++i)
				{
				alleles.add(this.alleleBinding.entryToObject(in));
				}
			GenotypeBuilder gb=new GenotypeBuilder(sampleName,alleles);
			gb.GQ(in.readInt());
			gb.DP(in.readInt());
			return gb.make();
			}
		@Override
		public void objectToEntry(Genotype g, TupleOutput w)
			{
			w.writeString(g.getSampleName());
			List<Allele> alleles= g.getAlleles();
			w.writeInt(alleles.size());
			for(int i=0;i< alleles.size();++i)
				{
				Allele a =alleles.get(i);
				this.alleleBinding.objectToEntry(a,w);
				}
			w.writeInt(g.hasGQ()?g.getGQ():-1);
			w.writeInt(g.hasDP()?g.getDP():-1);
			}
		}

	
	private static class ChromPosRefBinding
		extends TupleBinding<ChromPosRef>
		{
		AlleleBinding alleleBinding=new AlleleBinding();
		@Override
		public ChromPosRef entryToObject(TupleInput in)
			{
			ChromPosRef o=new ChromPosRef();
			o.chrom=in.readString();
			o.pos1=in.readInt();
			o.ref = alleleBinding.entryToObject(in);
			return o;
			}
		@Override
		public void objectToEntry(ChromPosRef o, TupleOutput out) {
			out.writeString(o.chrom);
			out.writeInt(o.pos1);
			alleleBinding.objectToEntry(o.ref, out);
			}
		}
	private final ChromPosRefBinding chromPosRefBindingInstance=new ChromPosRefBinding();
	
	
	
	/** class Allele
	
	/** coverage for position on genome */
	private class GenotypeContextBinding
		extends TupleBinding<GenotypesContext>
		{
		private GenotypeBinding genotypeBinding=new GenotypeBinding();
		@Override
		public GenotypesContext entryToObject(TupleInput in)
			{
			int n=in.readInt();
			ArrayList<Genotype> genotypes=new ArrayList<Genotype>(n);
			for(int i=0;i<n;++i) genotypes.add(this.genotypeBinding.entryToObject(in));
			return GenotypesContext.create(genotypes);
			}
		@Override
		public void objectToEntry(GenotypesContext o, TupleOutput w)
			{
			ArrayList<Genotype> genotypes=new ArrayList<Genotype>();
			for(int i=0;i< o.size();++i)
				{
				Genotype g=o.get(i);
				switch(g.getType())
					{
					case HET: case HOM_REF: case HOM_VAR: 
						{
						genotypes.add(g);
						break;
						}
					default:break;
					}
				}
			w.writeInt(genotypes.size());
			for(Genotype g:genotypes)
				{
				this.genotypeBinding.objectToEntry(g,w);
				}
			}
		}
	
	private  final GenotypeContextBinding genotypeContextBindingInstance=new GenotypeContextBinding();
	
	
	
	private long ID_GENERATOR=0L;
	
	
	
	private abstract class TreeNode
		{
		long nodeid=(++ID_GENERATOR);
		long score=0L;
		private Set<String> _cacheSamples=null;
		
		abstract void collectSampleIds(Set<String> set);
		public Set<String> getSamples()
			{
			if(_cacheSamples==null)
				{
				_cacheSamples=new HashSet<String>();
				collectSampleIds(_cacheSamples);
				}
			return _cacheSamples;
			}
		
		protected abstract Counter<GenotypeType> getGenotypeTypes(final GenotypesContext ctx);
		
		long distance(GenotypeType t1,GenotypeType t2)
			{
			if(t1.equals(t2)) return 0;
			
			if((t1==GenotypeType.HET &&
				(t2.equals(GenotypeType.HOM_REF) || t2.equals(GenotypeType.HOM_VAR)))
				)
					{
					return 5L;
					}
			if(t1==GenotypeType.HOM_REF && t2==GenotypeType.HOM_VAR)
				{
				return 15L;
				}
			return distance(t2,t1);
			}
		
		long distance(TreeNode other,final GenotypesContext ctx)
			{
			Counter<GenotypeType> x1= this.getGenotypeTypes(ctx);
			Counter<GenotypeType> x2= other.getGenotypeTypes(ctx);
			if(x1.isEmpty() && x2.isEmpty()) return 0L;
			if(x1.isEmpty() || x2.isEmpty()) return 1L;
			GenotypeType t1= x1.getMostFrequent();
			GenotypeType t2= x1.getMostFrequent();
			
			return distance(t1,t2);
			}

		
		protected GenotypeType getGenotypeType(Genotype g)
			{
			if(g==null ) return GenotypeType.UNAVAILABLE;
			GenotypeType gt=g.getType();		
			return ( gt.equals(GenotypeType.MIXED) ||
					gt.equals(GenotypeType.NO_CALL) ?
					GenotypeType.UNAVAILABLE:
					gt);
			}
		
		boolean containsAnySampleIds(Set<String> set)
			{
			for(String s: getSamples())
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
		
		}
	
	private class OneSample extends TreeNode
		{
		private String sample;
		OneSample(String sample)
			{
			this.sample = sample;
			}
		
		@Override
		protected Counter<GenotypeType> getGenotypeTypes(final GenotypesContext ctx)
			{
			Counter<GenotypeType> c=new Counter<GenotypeType>();
			GenotypeType t = getGenotypeType(ctx.get(this.sample));
			switch(t)
				{
				case HOM_REF: 
				case HOM_VAR:
				case HET: c.incr(t);break;
				default:break;
				}
			return c;
			}
		
		@Override
		void collectSampleIds(Set<String> set)
			{
			set.add(sample);
			}
		
		@Override
		void writeGraphizDot(PrintStream out) {
			out.println(getNodeId()+"[label=\""+this.sample+"\"];");
			}
		@Override
		 void writeGexfNodes(XMLStreamWriter w) throws XMLStreamException
			{
			w.writeEmptyElement("node");
				w.writeAttribute("id",getNodeId());
				w.writeAttribute("label",this.sample);
			}
		
		@Override
		void writeGexfEdges(XMLStreamWriter w) throws XMLStreamException
			{
			//empty
			}
		@Override
		void writeNewick(PrintStream out)
			{
			out.print(this.sample);
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
		void collectSampleIds(Set<String> set) {
			t1.collectSampleIds(set);
			t2.collectSampleIds(set);
			}
		
		
		@Override
		protected Counter<GenotypeType> getGenotypeTypes(final GenotypesContext ctx)
			{
			Counter<GenotypeType> c=new Counter<GenotypeType>();
			c.putAll(t1.getGenotypeTypes(ctx));
			c.putAll(t2.getGenotypeTypes(ctx));
			/* comparing two samples but two diferent genotypes
			 * keep the one with the highest QUAL
			 */
			if(c.getCountCategories()==2 && 
				t1 instanceof OneSample &&
				t2 instanceof OneSample )
				{
				Genotype g1= ctx.get(OneSample.class.cast(t1).sample);
				Genotype g2= ctx.get(OneSample.class.cast(t2).sample);
				if(g1.hasGQ() && g2.hasGQ())
					{
					c=new Counter<GenotypeType>();
					c.incr(g1.getGQ()>g2.getGQ()?
						getGenotypeType(g1):
						getGenotypeType(g2)
						);
					}
				}
			return c;
			}

		@Override
		void writeGraphizDot(PrintStream out)
			{
			t1.writeGraphizDot(out);
			t2.writeGraphizDot(out);
			out.println(this.getNodeId()+"[shape=point];");
			out.println(t1.getNodeId()+" -- "+ this.getNodeId()+";");
			out.println(t2.getNodeId()+" -- "+ this.getNodeId()+";");
			
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
			for(String nid:new String[]{
					this.t1.getNodeId(),
					this.t2.getNodeId()}
					)
				{
				w.writeEmptyElement("edge");
				w.writeAttribute("id","E"+this.nodeid+"_"+nid);
				w.writeAttribute("type","directed");
				w.writeAttribute("source",nid);
				w.writeAttribute("target",this.getNodeId());
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
		}
	

	
	/** map TidPos to CovRow */
	private Environment env;
	/** map TidPos to CovRow */
	private Database pos2cov;
	/** all samples */
	private Set<String> samples=new HashSet<String>();
	/** count snps */
	private long snp_count=0L;
	/** ignore variant if genotype missing */
	private boolean ignore_if_genotype_missing=false;
	
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
			GenotypesContext row = this.genotypeContextBindingInstance.entryToObject(data);
			
			if(this.ignore_if_genotype_missing)
				{
				Counter<GenotypeType> counter=new Counter<>();
				boolean ok=true;
				for(String sample:this.samples)
					{
					Genotype g=row.get(sample);
					if(g==null) { ok=false; break;}
					switch(g.getType())
						{
						case NO_CALL:
						case UNAVAILABLE: ok=false;break;
						default:break;
						}
					if(!ok) break;
					counter.incr(g.getType());
					}
				
				if(ok && counter.getCountCategories()<2)
					{
					ok=false;
					}
				
				if(!ok)
					{
					++count_deleted;
					cursor.delete();
					continue;
					}
				}
			}
		cursor.close();
		info("After cleanup, removed "+count_deleted+" snps.");
		}
	private TreeNode matrix(Transaction txn, List<TreeNode> parents)
		{
		long last=System.currentTimeMillis();
		
		/* faster computing ? , parse BerkeleyDB one time */
		if(Math.pow(parents.size(),2)/2.0 <500000.0)
			{
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
			Cursor cursor=null;
			cursor=this.pos2cov.openCursor(txn, null);
			long nVar=0L;
			while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				long now=System.currentTimeMillis();
				if(now-last > 10*1000)
					{
					info("iteration "+nVar+"/"+this.snp_count+" for "+parents.size());
					last=now;
					}
				
				GenotypesContext row = this.genotypeContextBindingInstance.entryToObject(data);
				for(MergedNodes merged : nodes)
					{
					merged.score += merged.t1.distance(merged.t2,row);
					}
				++nVar;
				}
			cursor.close();
			
			Collections.sort(nodes,new Comparator<MergedNodes>()
				{
				public int compare(MergedNodes o1, MergedNodes o2)
					{
					if(o1.score< o2.score) return -1;
					if(o1.score> o2.score) return 1;
					return 0;
					}
				});
			return nodes.get(0);
			}
		else
			{
			MergedNodes best=null;
			for(int x=0;x< parents.size();++x)
				{
				TreeNode tn1= parents.get(x);
				
				
				for(int y=x+1;y< parents.size();++y)
					{
					long now=System.currentTimeMillis();
					if(now-last > 10*1000)
						{
						info("iteration ["+(x+1)+","+(y+1)+"]/"+parents.size());
						last=now;
						}
					
					TreeNode tn2= parents.get(y);
					DatabaseEntry key=new DatabaseEntry();
					DatabaseEntry data=new DatabaseEntry();
					MergedNodes merged= new MergedNodes(tn1, tn2);
					Cursor cursor=null;
					cursor=this.pos2cov.openCursor(txn, null);
					while(cursor.getNext(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
						{
						GenotypesContext row = this.genotypeContextBindingInstance.entryToObject(data);
						
						merged.score += tn1.distance(tn2,row);
						}
					cursor.close();
					if(best==null || best.score > merged.score)
						{
						best=merged;
						}
					}
				}
			return best;
			}
		
		}
	
	private TreeNode matrix(Transaction txn)
		{
		ArrayList<TreeNode> nodes=new ArrayList<TreeNode>();
		//initialize population
		for(String sample:this.samples)
			{
			nodes.add(new OneSample(sample));
			}
		
		while(nodes.size()>1)
			{
			info("nodes.count= "+nodes.size());
			if(nodes.size()==2)
				{
				return new MergedNodes(
						nodes.get(0),
						nodes.get(1)
						);
				}
			
			
			TreeNode best = matrix(txn, nodes);
			Set<String> bestsamples= best.getSamples();
			
			ArrayList<TreeNode> newnodes=new ArrayList<TreeNode>();

			for(TreeNode tn:nodes)
				{
				if(tn.containsAnySampleIds(bestsamples)) continue;
				newnodes.add(tn);
				}
			newnodes.add(best);
			
			nodes=newnodes;
			}
		return nodes.get(0);
		}
	
	
	private void readvcf(Transaction txn,VcfIterator in)
		{
		Cursor cursor=null;
		VCFHeader header = in.getHeader();
		if(!header.hasGenotypingData()) return;
		this.samples.addAll(header.getSampleNamesInOrder());
		cursor= this.pos2cov.openCursor(txn, null);

		DatabaseEntry key = new DatabaseEntry();
		DatabaseEntry data = new DatabaseEntry();
		ChromPosRef chromPosRef=new ChromPosRef();
		while(in.hasNext())
			{
			VariantContext ctx = in.next();
			
			chromPosRef.chrom=ctx.getChr();
			chromPosRef.pos1=ctx.getStart();
			chromPosRef.ref=ctx.getReference();
			
			chromPosRefBindingInstance.objectToEntry(chromPosRef, key);
			
			if(cursor.getSearchKey(key, data, LockMode.DEFAULT)==OperationStatus.SUCCESS)
				{
				GenotypesContext row = genotypeContextBindingInstance.entryToObject(data);
				
				for(int i=0;i< ctx.getNSamples();++i)
					{
					Genotype g = ctx.getGenotype(i);
					if(!g.isCalled()) continue;
					Genotype old= row.get(g.getSampleName());
					if( old!=null && 
						old.hasGQ() && g.hasGQ() &&
						old.getGQ() < g.getGQ() 	
						)
						{
						continue;
						}
					
					row.add(g);
					}
				genotypeContextBindingInstance.objectToEntry(row, data);
				cursor.putCurrent(data);
				}
			else
				{
				genotypeContextBindingInstance.objectToEntry(ctx.getGenotypes(), data);
				cursor.put(key,data);
				}
			
			}
		cursor.close();
		}
	

	
	
	@Override
	public String getProgramDescription() {
		return "";
		}
	private enum OUT_FMT
		{
		dot,gexf,newick
		};
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -B bdb home. "+getMessageBundle("berkeley.db.home"));
		out.println(" -f (format) one of "+Arrays.toString(OUT_FMT.values()));
		//out.println(" -i ignore SNP if genotype missing"); TODO
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		OUT_FMT outformat=OUT_FMT.dot;
		File bdbHome=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"B:f:i"))!=-1)
			{
			switch(c)
				{
				case 'i': this.ignore_if_genotype_missing=true;break;
				case 'B': bdbHome=new File(opt.getOptArg());break;
				case 'f': outformat=OUT_FMT.valueOf(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(bdbHome==null)
			{
			error("undefined bdb home");
			return -1;
			}
		
		
		Transaction txn=null;
		VcfIterator iter=null;
		try
			{
			info("Opening "+bdbHome);
			EnvironmentConfig cfg=new EnvironmentConfig();
			cfg.setAllowCreate(true);
			cfg.setReadOnly(false);
			this.env= new Environment(bdbHome, cfg);
			
			DatabaseConfig cfg2=new DatabaseConfig();
			cfg2.setAllowCreate(true);
			cfg2.setTemporary(true);
			cfg2.setReadOnly(false);
			this.pos2cov=this.env.openDatabase(txn, "vcf", cfg2);
			
			
			if(opt.getOptInd()==args.length)
				{
				info("reading from stdin");
				iter=VCFUtils.createVcfIteratorStdin();
				readvcf(txn,iter);
				iter.close();
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("reading "+filename);
					iter=VCFUtils.createVcfIterator(filename);
					readvcf(txn,iter);
					iter.close();
					}
				}
			
			if(this.samples.size()<2)
				{
				error("Not enoug samples");
				return -1;
				}
			
			info("SNPs: "+this.snp_count);
			info("Samples: "+this.samples.size());
			
			
			cleanup(txn);
			
			this.snp_count= this.pos2cov.count();
			
			if(this.snp_count==0L)
				{
				error("Not enough SNPS.");
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
					  w.writeCharacters(getProgramName()+" by "+getAuthorName());
					  w.writeEndElement();
					 
					  w.writeStartElement("description");
					  w.writeCharacters(getProgramCommandLine());
					  w.writeEndElement();
					
					w.writeEndElement();//meta
					  
					  w.writeStartElement("graph");
					  w.writeAttribute("mode", "static");
					  w.writeAttribute("defaultedgetype", "directed");
					  
					  w.writeEmptyElement("attributes");
					  w.writeAttribute("class", "node");
					  w.writeAttribute("mode", "static");
					  
					  
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
		catch(Exception err)
			{
			error(err);
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
