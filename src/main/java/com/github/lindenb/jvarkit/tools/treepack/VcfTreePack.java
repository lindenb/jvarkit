package com.github.lindenb.jvarkit.tools.treepack;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class VcfTreePack extends  AbstractVcfTreePack
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfTreePack.class);

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
				this.get(ctx.getContig()).watch(ctx);
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
					@SuppressWarnings("unchecked")
					final BranchNode b=(BranchNode)parent;
					if(b instanceof SampleNodeFactory.MyNode)
						{
						final String sample= b.findKeyByValue(me);
						if(sample==null) throw new IllegalStateException();
						return sample;
						}
					me=parent;
					parent=parent.getParent();
					}
				throw new RuntimeException("Bad path: you MUST use a 'sample' before using a "+getName());
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
			 progress.watch(rec.getContig(),rec.getStart());
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
	public Collection<Throwable> call() throws Exception {
		if(super.listFactories) 
			{
			super.printAvailableFactories();
			return RETURN_OK;
			}
		setDimension(super.dimensionStr);
		buildFactoryChain(super.chainExpression);
		if(super.nodeFactoryChain.isEmpty())
			{
			return wrapException("no path defined");
			}
		
		VcfIterator in=null;
		final List<String> args= super.getInputFiles();
		try
			{
			if(args.isEmpty())
				{
				LOG.info("Reading stdin");
				in=VCFUtils.createVcfIteratorFromStream(stdin());
				scan(in);
				CloserUtil.close(in);
				}
			else
				{
				for(final String f: args)
					{
					in=VCFUtils.createVcfIterator(f);
					scan(in);
					CloserUtil.close(in);
					}
				}
			this.layout();
			this.svg();
			return RETURN_OK;
			}
		catch (final Exception err)
			{
			return wrapException(err);
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
