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
package com.github.lindenb.jvarkit.tools.burden;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.script.ScriptException;
import javax.script.SimpleBindings;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Node;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory.AFExtractor;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
/**
BEGIN_DOC

## The XML config

root is `<config>`

under config , we can find some `<filter>` , `<maf>` and `<handlers>`

A filter looks like `<filter id="abc">javascript expression filter the data </filter>` with the variable `variant`.
The script should return 'true' to accept the variant.


 `<maf id="abcd"/>` can have an @attribute `attribute` that will be the INFO field that contains the attribute
 otherwise it's considered a MAF extractor counting the genotypes. This attribute is parsed by AFExtractorFactory. So it can be 'AC/AN' or 'AF'

 `<handlers name="abcd"/>`  can have a `<filter>` or a reference to a filter `<filter ref='filter-id'>`
 `<handlers name="abcd"/>`  should have a `<case>` and a `<"ctrl">`
 
 
 `<case>` and a `<"ctrl">` both tags can have a @ref linking to a  `<maf>`. If there is no ref, then we use a MAF extractor counting the genotypes.


```xml

```

## Example



END_DOC
 */
@Program(
	name="casectrlplot",
	description="Plot CASE/CTRL data from VCF files",
	keywords={"maf","burden","case","control","plot","chart","vcf"},
	generate_doc=false
	)
public class CaseControlPlot extends Launcher
	{
	private static final Logger LOG=Logger.build(CaseControlPlot.class).make();
	@Parameter(names={"-o","--out"},description="Output Directory, or a filename ending with '.zip'",required=true)
	private File outputDirOrZip=null;
	@Parameter(names={"-tee","--tee"},description="Output the incoming VCF to stdout. Useful to capture intermediate results in pipelines.")
	private boolean teeToStdout=false;
	@Parameter(names={"-format","--format"},description="How to print doubles, printf-life precision format")
	private String precisionFormat="%.5f";
	@Parameter(names={"-prefix","--prefix"},description="Output files prefix")
	private String prefix="";
	@Parameter(names={"-ped","--pedigree"},description=Pedigree.OPT_DESCRIPTION+" If not defined, I will try to extract the pedigree from the VCFheader.")
	private File pedigreeFile;
	@Parameter(names={"-c","--config"},description="XML config file",required=true)
	private File xmlFile=null;
	
	
	private static interface MafExtractor
		{
		public Double apply(final VariantContext ctx,Allele alt,final Set<Pedigree.Person> persons);
		}
	
	private static class DefaultMafExtractor implements MafExtractor
		{
		private final AFExtractor afExtractor;
		public DefaultMafExtractor(final AFExtractor afExtractor) {this.afExtractor = afExtractor;}
		
		@Override
		public Double apply(final VariantContext ctx,final Allele alt,final Set<Pedigree.Person> persons) {
			final List<Double> L= this.afExtractor.parse(ctx);
			int index=ctx.getAlleleIndex(alt);
			if( index<=0) return null;//can't be REF==0
			if(index>=L.size()) return null;
			return L.get(index);
			}
		}
	
	private static class GenotypeMafExtractor implements MafExtractor
		{
		@Override
		public Double apply(final VariantContext ctx,final Allele alt,final Set<Pedigree.Person> samples) {
			final MafCalculator calc=new MafCalculator(alt, ctx.getContig());
			for(final Pedigree.Person sample: samples) {
				calc.add(ctx.getGenotype(sample.getId()), sample.isMale());
				}
			return calc.isEmpty()?null:calc.getMaf();
			}
		}
	
	
	
	/** interface to extract X/Y data from a variant */
	private class CaseControlExtractor
		{
		private String name;
		private Predicate<VariantContext> variantPredicate = C -> true;
		private MafExtractor caseExtractor = new GenotypeMafExtractor();
		private MafExtractor ctrlExtractor= new GenotypeMafExtractor();
		private PrintWriter pw;

		/** name of this extractor */
		public String getName() {
			return this.name;
			}
		/** Variant context filter: should I accept this variant ? */
		public Predicate<VariantContext> getFilter()
			{
			return this.variantPredicate;
			}
		/** case extractor */
		public MafExtractor getCaseMAFExtractor() {
			return this.caseExtractor;
			}
		/** control extractor */
		public MafExtractor getControlMAFExtractor() {
			return this.ctrlExtractor;
			}

		public  void visit(
				final VariantContext vc,
				final Set<Pedigree.Person> cases,
				final Set<Pedigree.Person> controls
				)
			{
			if(vc==null || !vc.isVariant()) return;
			if(getFilter()!=null)
				{
				if(!getFilter().test(vc)) return ;
				}
			final List<Allele> alts = vc.getAlternateAlleles();
			for(int altidx=0; altidx < alts.size();++altidx)
				{
				final Allele alt= alts.get(altidx);
				final Double casex = this.getCaseMAFExtractor().apply(vc, alt,cases);
				if(casex==null || casex<0.0 || casex>1.0) continue;
				final Double ctrly = this.getControlMAFExtractor().apply(vc, alt,controls);
				if(ctrly==null || ctrly<0.0 || ctrly>1.0) continue;
				
				pw.printf(precisionFormat,casex);
				pw.print('\t');
				pw.printf(precisionFormat,ctrly);
				pw.println();
				}
			}
		public void close() {
			if(pw!=null) pw.flush();
			CloserUtil.close(pw);
			}
		}
	
	private class JSPredigate implements Predicate<VariantContext>
		{
		final javax.script.CompiledScript compiled;
		final SimpleBindings bindings = new SimpleBindings();
		JSPredigate(final VCFHeader header,final String expression) throws Exception
			{
			this.compiled = compileJavascript(expression, null);
			this.bindings.put("header",header);
			this.bindings.put("tools",new VcfTools(header));
			}
		
		@Override
		public boolean test(final VariantContext t)
			{
			this.bindings.put("variant",t);
			Object o =null;
			try
				{
				o= this.compiled.eval(this.bindings);
				}
			catch (ScriptException e)
				{
				throw new RuntimeException(e);
				}
			if(o!=null && o instanceof Boolean)
				{
				return Boolean.class.cast(o);
				}
			throw new RuntimeException("Script returned somethng that is a not a boolean");
			}
		}	
	
	
	private List<CaseControlExtractor> parseConfigFile(final VCFHeader header) throws Exception {
		
		Document dom = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(this.xmlFile);
		Element root = dom.getDocumentElement();
		if(root==null || !root.getNodeName().equals("config"))
				throw new JvarkitException.XmlDomError(root,"Root is note <config>");
		final Map<String,JSPredigate> id2variantFilter = new HashMap<>();
		
		// first pass, collect filters
		for(Node n1=root.getFirstChild();n1!=null;n1=n1.getNextSibling())
			{
			if(n1.getNodeType()!=Node.ELEMENT_NODE) continue;
			Element e1=Element.class.cast(n1);
			if(!e1.getNodeName().equals("filter")) continue;
			final String filterid=e1.getAttribute("id");
			if(StringUtil.isBlank(filterid)) throw new JvarkitException.XmlDomError(e1,"@id missing");
			if(id2variantFilter.containsKey(filterid))  throw new JvarkitException.XmlDomError(e1,"duplicate filter id : "+filterid);
			final String expression = e1.getTextContent();
			if(StringUtil.isBlank(expression)) throw new JvarkitException.XmlDomError(e1,"expression missing");
			id2variantFilter.put(filterid,new JSPredigate(header, expression));
			}
		
		
		//second pass collect maf extractors
		final Map<String,MafExtractor> id2mafExtractor = new HashMap<>();
		for(Node n1=root.getFirstChild();n1!=null;n1=n1.getNextSibling())
			{
			if(n1.getNodeType()!=Node.ELEMENT_NODE) continue;
			Element e1=Element.class.cast(n1);
			if(!e1.getNodeName().equals("maf")) continue;
			final String mafid =e1.getAttribute("id");
			if(StringUtil.isBlank(mafid)) throw new JvarkitException.XmlDomError(e1,"@id missing");
			if(id2mafExtractor.containsKey(mafid))  throw new JvarkitException.XmlDomError(e1,"duplicate maf id : "+mafid);
			final MafExtractor mafExtractor; 
			if(e1.hasAttribute("attribute"))
				{
				final String tag = e1.getAttribute("attribute");
				if(StringUtil.isBlank(tag)) throw new JvarkitException.XmlDomError(e1,"@attribute is empty");
				final List<AFExtractor> extractors = new AFExtractorFactory().parseFieldExtractors(tag);
				if(extractors.size()!=1) {
					throw new JvarkitException.XmlDomError(e1,"expected one AF extractor but got "+extractors);
				}
				final DefaultMafExtractor at= new DefaultMafExtractor(extractors.get(0));
				mafExtractor = at;
				}
			else
				{
				final GenotypeMafExtractor at= new GenotypeMafExtractor();
				mafExtractor = at;
				}
			id2mafExtractor.put(mafid, mafExtractor);
			}
		final List<CaseControlExtractor> excractors= new ArrayList<>();
		//parse handlers
		for(Node n1=root.getFirstChild();n1!=null;n1=n1.getNextSibling())
			{
			if(n1.getNodeType()!=Node.ELEMENT_NODE) continue;
			Element e1=Element.class.cast(n1);
			if(!e1.getNodeName().equals("handler")) continue;
			final CaseControlExtractor extractor = new CaseControlExtractor();
			extractor.name  = e1.getAttribute("name");
			if(StringUtil.isBlank(extractor.name)) throw new JvarkitException.XmlDomError(e1,"@name missing");
			
			extractor.name  = this.prefix + extractor.name;
			
			if(excractors.stream().filter(F->F.getName().equals(extractor.name)).findAny().isPresent()) {
				throw new JvarkitException.XmlDomError(e1,"duplicate name"+extractor.name);
				}
			
			for(Node n2=n1.getFirstChild();n2!=null;n2=n2.getNextSibling())
				{
				if(n2.getNodeType()!=Node.ELEMENT_NODE) continue;
				Element e2=Element.class.cast(n2);
				if(e2.getNodeName().equals("filter")) {
					final Predicate<VariantContext> expr;
					if(e2.hasAttribute("ref"))
						{
						final String filterid = e2.getAttribute("ref");
						JSPredigate predicate = id2variantFilter.get(filterid); 
						
						if(predicate==null)
							{
							throw new JvarkitException.XmlDomError(e2,"no such filter:"+filterid);
							}
						expr =predicate;
						}
					else
						{
						final String expressionstr = e2.getTextContent();
						if(StringUtil.isBlank(expressionstr)) throw new JvarkitException.XmlDomError(e2,"expression missing");
						expr = new JSPredigate(header,expressionstr);
						}
					extractor.variantPredicate = extractor.variantPredicate.and(expr);
					}
				else if(e2.getNodeName().equals("case") || 
						e2.getNodeName().equals("ctrl"))
					{
					final MafExtractor mafextractor;
					if(e2.hasAttribute("ref"))
						{
						final String mafid = e2.getAttribute("ref");
						if(!id2mafExtractor.containsKey(e2.getAttribute(mafid)))
							{
							throw new JvarkitException.XmlDomError(e2,"no such mafextractor:"+mafid);
							}
						mafextractor = id2mafExtractor.get(mafid);
						}
					else
						{
						final GenotypeMafExtractor genotypeMafExtractor=new GenotypeMafExtractor();
						mafextractor =genotypeMafExtractor;
						}
					if(e2.getNodeName().equals("case"))
						 {
						 extractor.caseExtractor = mafextractor;
						 }
					else 
						 {
						 extractor.ctrlExtractor = mafextractor;
						 }
					}
				else
					{
					LOG.error("unknown XML element "+e2.getNodeName());
					}
				}
			excractors.add(extractor);
			}
		return excractors;
		}

	@Override
	public int doWork(final List<String> args)
		{
		ArchiveFactory archiveFactory=null;
		VCFIterator in = null;
		VariantContextWriter teeVariantWriter = null;
		final List<CaseControlExtractor> excractors = new ArrayList<>();
		try
			{
			
			in = super.openVCFIterator(oneFileOrNull(args));
			final VCFHeader header= in.getHeader();
			excractors.addAll( parseConfigFile(header) ) ;
			final Pedigree pedigree;
			if( this.pedigreeFile!=null) {
				pedigree = Pedigree.newParser().parse(this.pedigreeFile);
				}
			else
				{
				pedigree = Pedigree.newParser().parse(header);
				}
			if(pedigree==null || pedigree.isEmpty()) {
				LOG.error("No pedigree defined , or it is empty");
				return -1;
				}
			final Set<Pedigree.Person> casepersons = pedigree.getPersons().
					stream().
					filter(F->F.isAffected() && header.getSampleNamesInOrder().indexOf(F.getId())!=-1).
					collect(Collectors.toSet());
					
			if(casepersons.isEmpty()){
					LOG.error("No Affected individuals in pedigree/header");
					return -1;
					}
			
			final Set<Pedigree.Person> controlpersons = pedigree.getPersons().
					stream().
					filter(F->F.isUnaffected() && header.getSampleNamesInOrder().indexOf(F.getId())!=-1).
					collect(Collectors.toSet());
			
			if(controlpersons.isEmpty()){
					LOG.error("No Unaffected individuals in pedigree/header");
					return -1;
					}
			if( this.teeToStdout)
				{
				teeVariantWriter = super.openVariantContextWriter(null);
				teeVariantWriter.writeHeader(header);
				}
			archiveFactory = ArchiveFactory.open(this.outputDirOrZip);
			for(final CaseControlExtractor extractor : excractors) {
				extractor.pw = archiveFactory.openWriter(extractor.name);
				}
			final SAMSequenceDictionaryProgress progress= new SAMSequenceDictionaryProgress(header);
			while(in.hasNext())
				{
				final VariantContext ctx = progress.watch(in.next());
				if(teeVariantWriter!=null) teeVariantWriter.add(ctx);
				for(final CaseControlExtractor handler : excractors) {
					handler.visit(ctx,casepersons,controlpersons);
					}
				}
			for(final CaseControlExtractor handler : excractors) {
				handler.close();
				}
			progress.finish();
			if(teeVariantWriter!=null) {
				teeVariantWriter.close();
				teeVariantWriter=null;
				}
			in.close();in=null;
			archiveFactory.close();archiveFactory=null;
			return RETURN_OK;
			}
		catch (final Exception e)
			{
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(archiveFactory);
			CloserUtil.close(teeVariantWriter);
			for(final CaseControlExtractor handler : excractors) {
				handler.close();
				}
			}
		}
	
	public static void main(final String[] args)
		{
		new CaseControlPlot().instanceMainWithExit(args);
		}

	}
