package com.github.lindenb.jvarkit.tools.treepack;

import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfTreePack extends AbstractTreePackCommandLine<VariantContext>
	{
	/** Group by Chromosome */
	private class ChromosomeNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(ChromosomeNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(VariantContext ctx)
				{
				this.get(ctx.getChr()).watch(ctx);
				}
			}
		@Override
		public String getDescription()
			{
			return "chromosomes";
			}
		@Override
		public String getName()
			{
			return "chrom";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}

	/** Group by Sample */
	private class SampleNodeFactory extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(SampleNodeFactory.this,parent,label);
				}
			
			@Override
			public void watch(VariantContext ctx)
				{
				if(!ctx.hasGenotypes()) return;
				for(Genotype g:ctx.getGenotypes())
					{
					if(!g.isAvailable()) continue;
					if(g.isNoCall()) continue;
					this.get(g.getSampleName()).watch(ctx);
					}
				
				}
			}
		@Override
		public String getDescription()
			{
			return "sample-name";
			}
		@Override
		public String getName()
			{
			return "sample";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}
	/** Group by Genotype */
	private abstract class AbstractSubSampleNodeFactory extends NodeFactory
		{
		private abstract class MyGNode extends BranchNode
			{
			private MyGNode(BranchNode parent,String label)
				{
				super(AbstractSubSampleNodeFactory.this,parent,label);
				}
			private String findSample()
				{
				AbstractNode me=this;
				AbstractNode parent=me.getParent();
				while(parent!=null)
					{
					if(parent.isLeaf()) break;
					BranchNode b=(BranchNode)parent;
					if(b instanceof SampleNodeFactory.MyNode)
						{
						String sample= b.findKeyByValue(me);
						if(sample==null) throw new IllegalStateException();
						return sample;
						}
					me=parent;
					parent=parent.getParent();
					}
				throw new PicardException("Bad path: you MUST use a 'sample' before using a "+getName());
				}
			public abstract void watch(Genotype g,VariantContext ctx);
			@Override
			public void watch(VariantContext ctx)
				{
				if(!ctx.hasGenotypes()) return;
				String sample=this.findSample();
				Genotype genotype=ctx.getGenotype(sample);
				if(genotype==null) throw new IllegalStateException();
				watch(genotype,ctx);
				}
			}
		}
	
	
	/** Group by Genotype */
	private class GenotypeNodeFactory extends AbstractSubSampleNodeFactory
		{
		private class MyNode extends AbstractSubSampleNodeFactory.MyGNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(parent,label);
				}
			
			@Override
			public void watch(Genotype g,VariantContext ctx)
				{
				this.get(g.getType().name()).watch(ctx);
				}
			}
		@Override
		public String getDescription()
			{
			return "genotype";
			}
		@Override
		public String getName()
			{
			return "genotype";
			}
		
		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}

	
	/** Group by MAPQ */
	private class QualityFactory
		extends NodeFactory
		{
		private class MyNode extends BranchNode
			{
			private MyNode(BranchNode parent,String label)
				{
				super(QualityFactory.this,parent,label);
				}
			
			@Override
			public void watch(VariantContext ctx)
				{
				String qualStr="*";
				if(ctx.hasLog10PError())
					{
					int qual=(int)ctx.getPhredScaledQual();
					qual=((int)(qual/10.0));
					qualStr=(int)(qual)*10+"-"+((qual+1)*10);
					}
				this.get(qualStr).watch(ctx);
				}
			}
		@Override
		public String getDescription()
			{
			return "quality (window.size=10)";
			}
		@Override
		public String getName() {
			return "qual";
			}

		@Override
		public BranchNode createBranch(BranchNode parent,String label)
			{
			return new MyNode(parent,label);
			}
		}

	
	 private void scan(VcfIterator iter)
		 {
		 VCFHeader header=iter.getHeader();
		 SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
		 while(iter.hasNext())
			 {
			 VariantContext rec=iter.next();
			 progress.watch(rec.getChr(),rec.getStart());
			 root.watch(rec);
			 }
		 progress.finish();
		 }

	private  List<NodeFactory> _factories=null;

	@Override
	protected List<NodeFactory> getAllAvailableFactories()
		{
		if(_factories==null)
			{
			_factories=new ArrayList<NodeFactory>();
			_factories.add(new ChromosomeNodeFactory());
			_factories.add(new QualityFactory());
			_factories.add(new SampleNodeFactory());
			_factories.add(new GenotypeNodeFactory());
			}
		return _factories;
		}
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/VcfTreePack";
		}
	@Override
	public String getProgramDescription()
		{
		return "Create a TreeMap from one or more VCF. Ouput is a SVG file.";
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()))!=-1)
			{
			switch(c)
				{
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(super.nodeFactoryChain.isEmpty())
			{
			error("no path defined");
			return -1;
			}
		
		VcfIterator in=null;
		
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading stdin");
				in=VCFUtils.createVcfIteratorStdin();
				scan(in);
				CloserUtil.close(in);
				}
			else
				{
				for(int optind =opt.getOptInd(); optind < args.length;++optind)
					{
					String f = args[optind];
					in=VCFUtils.createVcfIterator(f);
					scan(in);
					CloserUtil.close(in);
					}
				}
			this.layout();
			this.svg();
			return 0;
			}
		catch (Exception e)
			{
			error(e);
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
	public static void main(String[] args)
		{
		new VcfTreePack().instanceMainWithExit(args);
		}

	}
