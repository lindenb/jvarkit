package com.github.lindenb.jvarkit.util.vcf;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import javax.script.CompiledScript;
import javax.script.ScriptException;
import javax.script.SimpleBindings;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.Text;

import com.github.lindenb.jvarkit.annotproc.IncludeSourceInJar;
import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;

import htsjdk.variant.vcf.VCFHeader;

@IncludeSourceInJar
public class XmlFilterVariantFactory
	{
	protected final Map<String,NamedFilterImpl> id2variantFilter = new HashMap<>();

	private class RefFilter implements Predicate<VariantContext>
		{
		final String ref;
		RefFilter(final String ref)
			{
			this.ref= ref;
			}
		@Override
		public boolean test(final VariantContext t)
			{
			final NamedFilterImpl delegate  = id2variantFilter.get(this.ref);
			if(delegate==null) throw new RuntimeException("Cannot find filter id="+this.ref);
			return delegate.test(t);
			}
		}
	
	public static interface NamedFilter extends Predicate<VariantContext>
		{
		public String getName();
		}
	private static class NamedFilterImpl implements NamedFilter
		{
		final String name;
		final Predicate<VariantContext> filter;
		NamedFilterImpl(final String name,Predicate<VariantContext> filter) {
			this.name = name;
			this.filter=filter;
			}
		@Override
		public String getName()
			{
			return this.name;
			}
		@Override
		public boolean test(VariantContext t)
			{
			return filter.test(t);
			}
		}

	
	private final VCFHeader header;
	private final VcfTools vcfTools;
	public XmlFilterVariantFactory(final VCFHeader header) {
		this.header = header;
		this.vcfTools = new VcfTools(header);
		}
	protected List<Element> childElements(final Node root)
		{
		if(root==null) return Collections.emptyList();
		final List<Element> L = new ArrayList<>();
		for(Node n1=root.getFirstChild();n1!=null;n1=n1.getNextSibling())
			{
			if(n1.getNodeType()==Node.TEXT_NODE)
				{
				final String content = Text.class.cast(n1).getData();
				if( content==null || StringUtil.isBlank(content)) continue;
				throw new JvarkitException.XmlDomError(root,"node contains text");
				}
			if(n1.getNodeType()!=Node.ELEMENT_NODE) continue;
			L.add(Element.class.cast(n1));
			}
		return L;
		}	
	
	
	public void parse(final Document dom) {
		Element root = dom.getDocumentElement();
		if(root==null) throw new JvarkitException.XmlDomError(dom,"No root");
		// first pass, collect filters
		for(final Element e1: childElements(root))
			{
			final Predicate<VariantContext> p= _filter(e1);
			if( p == null ) continue;
			final String filterid=e1.getAttribute("id");
			if(!filterid.matches("[a-zA-Z_][_a-zA-Z0-9]+"))
				{
				 throw new JvarkitException.XmlDomError(e1,"Bad id , no space or strange things please "+filterid);
				}
			if(id2variantFilter.containsKey(filterid))  throw new JvarkitException.XmlDomError(e1,"duplicate filter id : "+filterid);
			id2variantFilter.put(filterid,new NamedFilterImpl(filterid, p));
			}
		}
	
	private Predicate<VariantContext> _filter(Element e1) {
		if(e1==null) {
			return null;
			}
		else if(e1.getNodeName().equals("js") || e1.getNodeName().equals("javascript")) {
			return _jsfilter(e1);
			}
		else if(e1.getNodeName().equals("and")) {
			return _and(e1);
			}
		else if(e1.getNodeName().equals("or")) {
			return _or(e1);
			}
		else if(e1.getNodeName().equals("not") || e1.getNodeName().equals("negate")) {
			return _negate(e1);
			}
		else if(e1.getNodeName().equals("filter")) {
			if(!e1.hasAttribute("ref")) throw new JvarkitException.XmlDomError(e1, "@ref missing");
			return new RefFilter(e1.getAttribute("ref"));
			}
		return null;
		}
	
	private Predicate<VariantContext> _jsfilter(final Element root) {
		final String expression = root.getTextContent();
		if(StringUtil.isBlank(expression)) throw new JvarkitException.XmlDomError(root,"expression missing");
		return compileJSPredigate(expression);
		}
	
	
	private Predicate<VariantContext> _and(final Element root) {
		Predicate<VariantContext> P = VC->true;
		boolean found=false;
		for(final Element e1: childElements(root)) {
			Predicate<VariantContext> p = _filter(e1);
			if(p==null) throw new JvarkitException.XmlDomError(e1,"bad element");
			P = P.and(p);
			found=true;
			}
		if(!found) throw new JvarkitException.XmlDomError(root,"nothing found");
		return P;
		}
	private Predicate<VariantContext> _or(final Element root) {
		Predicate<VariantContext> P = VC->false;
		boolean found=false;
		for(final Element e1: childElements(root)) {
			Predicate<VariantContext> p = _filter(e1);
			if(p==null) throw new JvarkitException.XmlDomError(e1,"bad element");
			P = P.or(p);
			}
		if(!found) throw new JvarkitException.XmlDomError(root,"nothing found");
		return P;
		}
	private Predicate<VariantContext> _negate(final Element root) {
		Predicate<VariantContext> P = null;
		for(final Element e1: childElements(root)) {
			Predicate<VariantContext> pff = _filter(e1);
			if(pff==null) throw new JvarkitException.XmlDomError(e1,"bad element");
			if(P!=null) throw new JvarkitException.XmlDomError(root,"only one element allowed");
			P = pff.negate();
			}
		if(P==null) throw new JvarkitException.XmlDomError(root,"nothing found");
		return P;
		}

	private JSPredigate compileJSPredigate(final String jsExpr) {
		final javax.script.ScriptEngineManager manager = new javax.script.ScriptEngineManager();
		final javax.script.ScriptEngine engine = manager.getEngineByName("js");
		if(engine==null)
			{
			throw new JvarkitException.JavaScriptEngineNotFound();
			}
		final javax.script.Compilable compilingEngine = (javax.script.Compilable)engine;
		try {
			CompiledScript compiled = compilingEngine.compile(jsExpr);
			return new JSPredigate(compiled);
			}
		catch(final Exception err)
			{
			throw new RuntimeException(err);
			}
		}
	/** create simple bindings for javascript, insert 'header' and 'tools' 
	 * subclass have a chance to insert new things in the binding */
	public SimpleBindings createJavascriptBindings()  {
		final SimpleBindings bindings = new SimpleBindings();
		bindings.put("header",this.header);
		bindings.put("tools",this.vcfTools);
		return bindings;
		}
	
	private class JSPredigate implements Predicate<VariantContext>
		{
		final CompiledScript compiled;
		final SimpleBindings bindings = createJavascriptBindings();
		JSPredigate(final CompiledScript compiled) throws Exception
			{
			this.compiled = compiled;
			}
		
		@Override
		public boolean test(final VariantContext ctx)
			{
			this.bindings.put("variant",ctx);
			Object o =null;
			try
				{
				o= this.compiled.eval(this.bindings);
				}
			catch (final ScriptException e)
				{
				throw new JvarkitException.ScriptingError(e);
				}
			if(o!=null && o instanceof Boolean)
				{
				return Boolean.class.cast(o);
				}
			throw new JvarkitException.ScriptingError("Script returned somethng that is a not a boolean");
			}
		}	
	
	}
