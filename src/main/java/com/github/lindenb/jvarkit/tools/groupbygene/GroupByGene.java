package com.github.lindenb.jvarkit.tools.groupbygene;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFHeader;

import net.sf.picard.PicardException;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;


public class GroupByGene extends AbstractCommandLineProgram
	{
	private Set<String> sampleNames=new TreeSet<String>();
	private Set<String> user_gene_tags=new HashSet<String>();
	private SortingCollection<Call> sortingCollection=null;
	private boolean xml_output=false;
	
	private static class GeneName
		{
		String name;
		String type;
		GeneName(String name,String type)
			{
			this.name=name;
			this.type=type;
			}
		@Override
		public int hashCode()
			{
			final int prime = 31;
			int result = 1;
			result = prime * result +  name.hashCode();
			result = prime * result +  type.hashCode();
			return result;
			}
		@Override
		public boolean equals(Object o)
			{
			if (this == o) return true;
			if (o == null) return false;
			if (getClass() != o.getClass()) return false;
			GeneName g=(GeneName)o;
			return name.equals(g.name) && type.equals(g.type);
			}
		@Override
		public String toString() {
			return  name+"("+type+")";
			}
		
		}
	
	private static class Call
		{
		String chrom;
		int pos;
		String ref;
		GeneName gene;
		String sample;
		String posId;
		String a1=null;
		String a2=null;
		}
	
	private static class CallCodec
		extends AbstractDataCodec<Call>
		{
		@Override
		public void encode(DataOutputStream dos, Call c)
				throws IOException
			{
			dos.writeUTF(c.chrom);
			dos.writeInt(c.pos);
			dos.writeUTF(c.ref);
			dos.writeUTF(c.gene.name);
			dos.writeUTF(c.gene.type);
			
			dos.writeUTF(c.sample);
			dos.writeUTF(c.posId);
			dos.writeUTF(c.a1);
			dos.writeUTF(c.a2);
			}
		
		@Override
		public Call decode(DataInputStream dis) throws IOException
			{
			Call c= new Call();
			try {
				c.chrom=dis.readUTF();
			} catch (Exception e) {
				return null;
				}
			c.pos=dis.readInt();
			c.ref=dis.readUTF();
			String gName=dis.readUTF();
			String gType=dis.readUTF();
			c.gene=new GeneName(gName, gType);
			c.sample=dis.readUTF();
			c.posId=dis.readUTF();
			c.a1=dis.readUTF();
			c.a2=dis.readUTF();
			return c;
			}
		@Override
		public CallCodec clone() {
			return new CallCodec();
			}
		}
	
	private static class CallCmp
		implements Comparator<Call>
		{
		@Override
		public int compare(final Call o1, final Call o2)
			{
			int i=o1.chrom.compareTo(o2.chrom);
			if(i!=0) return i;
			 i= o1.gene.name.compareTo(o2.gene.name);
			 if(i!=0) return i;
			 i= o1.gene.type.compareTo(o2.gene.type);
			 return i;
			}
		}
	
	private SnpEffPredictionParser snpEffPredictionParser=null;
	private VepPredictionParser vepPredictionParser=null;
	
	private Set<GeneName> getGenes(VariantContext ctx)
		{
		HashSet<GeneName> set=new HashSet<GeneName>();
		for(VepPredictionParser.VepPrediction pred: this.vepPredictionParser.getPredictions(ctx))
			{
			String s=pred.getGeneName();
			if(s!=null)  set.add(new GeneName(s,"vep-gene-name"));
			s=pred.getEnsemblGene();
			if(s!=null)  set.add(new GeneName(s,"vep-ensembl-gene-name"));
			s=pred.getEnsemblTranscript();
			if(s!=null)
				{
				set.add(new GeneName(s,"vep-ensembl-transcript-name"));
				String ex=pred.getExon();
				if(ex!=null)
					{
					set.add(new GeneName(s+"|"+ex,"vep-ensembl-transcript-exon-name"));
					}
				}
			}
		for(SnpEffPredictionParser.SnpEffPrediction pred: this.snpEffPredictionParser.getPredictions(ctx))
			{
			String s=pred.getGeneName();
			if(s!=null)  set.add(new GeneName(s,"snpeff-gene-name"));
			s=pred.getEnsemblGene();
			if(s!=null)  set.add(new GeneName(s,"snpeff-ensembl-gene-name"));
			s=pred.getEnsemblTranscript();
			if(s!=null)  set.add(new GeneName(s,"snpeff-ensembl-transcript-name"));
			}
		for(String user_gene_tag:user_gene_tags)
			{
			if(user_gene_tag.isEmpty()) continue;
			if(user_gene_tag.equals(".")) continue;
			Object o=ctx.getAttribute(user_gene_tag);
			Set<String> tag=new HashSet<String>();
			if(o==null) 
				{
				//
				}
			else if(o.getClass().isArray())
				{
				for(Object o2:(Object[])o) tag.add(String.valueOf(o2));
				}
			else if(o instanceof java.util.Collection)
				{
				for(Object o2:(java.util.Collection<?>)o) tag.add(String.valueOf(o2));
				}
			else
				{
				tag.add(o.toString());
				}
			for(String t:tag)
				{
				if(t.isEmpty()) continue;
				set.add(new GeneName(t,"user:"+user_gene_tag));
				}
			}
		return set;
		}
	
	private void dump() throws IOException,XMLStreamException
		{
		XMLStreamWriter w=null;
		if(xml_output)
			{
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			w=xof.createXMLStreamWriter(System.out, "UTF-8");
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("genes");
			w.writeComment("Cmd line: "+this.getProgramCommandLine());
			w.writeComment("Version "+getVersion());
			w.writeStartElement("samples");
			w.writeAttribute("count", String.valueOf(sampleNames.size()));
			for(String sample:this.sampleNames)
				{
				w.writeStartElement("sample");
				w.writeCharacters(sample);
				w.writeEndElement();
				}
			w.writeEndElement();
			w.writeCharacters("\n");
			}
		else
			{
			System.out.print("#chrom");
			System.out.print('\t');
			System.out.print("min.POS");
			System.out.print('\t');
			System.out.print("max.POS");
			System.out.print('\t');
			System.out.print("gene.name");
			System.out.print('\t');
			System.out.print("gene.type");
			System.out.print('\t');
			System.out.print("samples.affected");
			System.out.print('\t');
			System.out.print("count.variations");
			for(String sample:this.sampleNames)
				{
				System.out.print('\t');
				System.out.print(sample);
				}
			System.out.println();
			}
		
		final CallCmp cmp=new CallCmp();
		List<Call> row=new ArrayList<Call>();
		CloseableIterator<Call> iter=sortingCollection.iterator();
		for(;;)
			{
			Call curr=null;
			if(iter.hasNext()) curr=iter.next();
			if(curr==null || (!row.isEmpty() && cmp.compare(curr, row.get(0))!=0))
				{
				if(!row.isEmpty())
					{
					int minPos=Integer.MAX_VALUE;
					int maxPos=Integer.MIN_VALUE;
					Set<String> affected=new HashSet<String>();
					Set<String> distinctMut=new HashSet<String>();
					Counter<String> sample2count=new Counter<String>();
					for(Call c:row)
						{
						minPos=Math.min(minPos, c.pos);
						maxPos=Math.max(maxPos, c.pos);
						sample2count.incr(c.sample);
						affected.add(c.sample);
						distinctMut.add(c.chrom+":"+c.pos+":"+c.ref);
						}
					Call first=row.get(0);
					if(w!=null)
						{
						w.writeStartElement("gene");
						w.writeAttribute("name", first.gene.name);
						w.writeAttribute("type", first.gene.type);
						
						w.writeAttribute("chrom", first.chrom);
						w.writeAttribute("min.POS",String.valueOf(minPos));
						w.writeAttribute("max.POS",String.valueOf(maxPos));
						w.writeAttribute("affected",String.valueOf(affected.size()));
						w.writeAttribute("variations",String.valueOf(distinctMut.size()));
						
						for(String sample:this.sampleNames)
							{
							if(sample2count.count(sample)==0L) continue;
							w.writeStartElement("sample");
							w.writeAttribute("name",sample);
							
							w.writeAttribute("count",String.valueOf(sample2count.count(sample)));
							for(Call c:row)
								{
								if(!c.sample.equals(sample)) continue;
								w.writeEmptyElement("genotype");
								w.writeAttribute("pos", String.valueOf(c.pos));
								w.writeAttribute("ref", c.ref);
								w.writeAttribute("A1", c.a1);
								w.writeAttribute("A2", c.a2);
								}
							w.writeEndElement();							
							}
						
						
						w.writeEndElement();//gene
						w.writeCharacters("\n");
						}
					else
						{
						System.out.print(first.chrom);
						System.out.print('\t');
						System.out.print(minPos);
						System.out.print('\t');
						System.out.print(maxPos);
						System.out.print('\t');
						System.out.print(first.gene.name);
						System.out.print('\t');
						System.out.print(first.gene.type);
						System.out.print('\t');
						System.out.print(affected.size());
						System.out.print('\t');
						System.out.print(distinctMut.size());
						for(String sample:this.sampleNames)
							{
							System.out.print('\t');
							System.out.print(sample2count.count(sample));
							}
						System.out.println();
						if(System.out.checkError()) break;
						}
					}
				if(curr==null) break;
				row.clear();
				}
			row.add(curr);
			}
		iter.close();
		if(w!=null)
			{
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			w.close();
			}
		}
	
	private void read(InputStream in)
		{
		VcfIterator iter=new VcfIterator(in);
		VCFHeader header=(VCFHeader)iter.getHeader();
		if(header.getSampleNamesInOrder()!=null)
			{
			this.sampleNames.addAll(header.getSampleNamesInOrder());
			}
		snpEffPredictionParser=new SnpEffPredictionParser(header);
		vepPredictionParser=new VepPredictionParser(header);
		SAMSequenceDictionary dict=header.getSequenceDictionary();
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict);
		while(iter.hasNext())
			{
			VariantContext ctx=iter.next();
			progress.watch(ctx.getChr(),ctx.getStart());
			for(GeneName g:getGenes(ctx))
				{
				for(Genotype genotype:ctx.getGenotypes())
					{
					if(!genotype.isAvailable()) continue;
					if(!genotype.isCalled()) continue;
					if(genotype.isNoCall()) continue;
					if(genotype.isHomRef()) continue;
					List<Allele> L=genotype.getAlleles();
					if(L==null || L.isEmpty()) continue;
					Call c=new Call();
					c.chrom=ctx.getChr();
					c.pos=ctx.getStart();
					c.ref=ctx.getReference().getDisplayString();
					c.gene=g;
					c.sample=genotype.getSampleName();
					c.posId=(ctx.getChr()+":"+ctx.getStart()+":"+ctx.getReference().getDisplayString()).toLowerCase();
					if(L.size()==1)
						{
						c.a1=genotype.getAllele(0).getDisplayString().toUpperCase();
						c.a2=c.a1;
						}
					else if(L.size()==2)
						{
						c.a1=genotype.getAllele(0).getDisplayString().toUpperCase();
						c.a2=genotype.getAllele(1).getDisplayString().toUpperCase();
						if(c.a1.compareTo(c.a2)>0)
							{
							String tmp=c.a1;
							c.a1=c.a2;
							c.a2=tmp;
							}
						}
					else
						{
						iter.close();
						throw new PicardException("cannot handle multi-ploidy "+ctx);
						}
					this.sortingCollection.add(c);
					}
				}	
			}
		iter.close();
		}
	
	
	

	@Override
	public String getProgramDescription() {
		return "Group VCF data by gene/transcript. By default it uses data from VEP , SnpEff";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -X XML output");
		out.println(" -T (tag) add Tag in INFO field containing the name of the genes.");
		
		super.printOptions(out);
		}
	
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/GroupByGene";
		}
	
	@Override
	public int doWork(String[] args)
		{
		int maxRecordsInRAM=10000;

		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"T:X"))!=-1)
			{
			switch(c)
				{
				case 'X': xml_output=true;break;
				case 'T': user_gene_tags.add(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		
		try
			{
			sortingCollection=SortingCollection.newInstance(
					Call.class,
					new CallCodec(),
					new CallCmp(),
					maxRecordsInRAM
					);
			sortingCollection.setDestructiveIteration(true);
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				read(System.in);
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					InputStream in=IOUtils.openURIForReading(filename);
					read(in);
					in.close();
					}
				}
			sortingCollection.doneAdding();

			info("Done reading. Now printing results.");

			dump();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			if(sortingCollection!=null) sortingCollection.cleanup();
			}
		}
	public static void main(String[] args) {
		new GroupByGene().instanceMainWithExit(args);
	}
}
