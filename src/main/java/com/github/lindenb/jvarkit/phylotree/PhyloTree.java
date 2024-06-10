/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.phylotree;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.variant.variantcontext.Allele;

public class PhyloTree {
public static final String DEFAULT_URI="https://raw.githubusercontent.com/genepi/phylotree-rcrs-17/main/src/tree.xml";
public static final String OPT_DESC="Path/URL to a valid XML formatted phylotree ( https://github.com/genepi/phylotree-rcrs-17 ).";
private final Map<String,HaploGroupImpl> name2haplogroup = new HashMap<>();

static public interface Polymorphism extends Comparable<Polymorphism> {
	public int getPosition();
	public Allele getAltAllele();
	}

private static class PolymorphismImpl implements Polymorphism {
	final int pos;
	final Allele alt;
	PolymorphismImpl(int pos, Allele alt) {
		this.pos = pos;
		this.alt = alt;
		}
	@Override
	public int getPosition() {
		return this.pos;
		}
	@Override
	public Allele getAltAllele() {
		return this.alt;
		}
	@Override
	public int compareTo(Polymorphism o) {
		int i = Integer.compare(getPosition(), o.getPosition());
		if(i!=0) return i;
		return getAltAllele().compareTo(o.getAltAllele());
		}
	@Override
	public int hashCode() {
		return Integer.hashCode(this.pos)%31 + this.alt.hashCode();
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof PolymorphismImpl)) return false;
		return this.pos==PolymorphismImpl.class.cast(obj).pos && 
			this.alt.equals(PolymorphismImpl.class.cast(obj).alt);
		}
	@Override
	public String toString() {
		return String.valueOf(pos)+this.alt.getDisplayString();
		}

	}

static public interface HaploGroup {
	public String getName();
	/** get polymporphisms including children */
	//public Set<Polymorphism> getAllPolymorphisms();
	public boolean hasDescendant(HaploGroup selfOrDescendant);
	}
static private class HaploGroupImpl implements HaploGroup {
	private final String sn;
	HaploGroupImpl parent=null;
	private final List<HaploGroupImpl> children = new ArrayList<>();
	private final Set<Polymorphism> polymorphisms = new HashSet<>();
	HaploGroupImpl(final String sn) {
		this.sn = sn;
		}
	private Set<Polymorphism> collectPolymorphisms(Set<Polymorphism> set) {
		set.addAll(this.polymorphisms);
		for(HaploGroupImpl hg:children) {
			hg.collectPolymorphisms(set);
			}
		return set;
		}
	
	@Override
	public boolean hasDescendant(HaploGroup selfOrDescendant) {
		if(this.equals(selfOrDescendant)) return true;
		for(HaploGroupImpl c: this.children) {
			if(c.hasDescendant(selfOrDescendant)) return true;
			}
	
		return false;
		}
	
	/* @Override
	public Set<Polymorphism> getAllPolymorphisms() {
		Set<Polymorphism> set = new TreeSet<>();
		return collectPolymorphisms(set);
		}*/
	@Override
	public String getName() {
		return sn;
		}
	@Override
	public int hashCode() {
		return sn.hashCode();
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof HaploGroup)) return false;
		return getName().equals(HaploGroup.class.cast(obj).getName());
		}
	@Override
	public String toString() {
		return sn;
		}
	}

/** get all haplogroups */
public Collection<HaploGroup> getAllHaploGroups() {
	return Collections.unmodifiableCollection(this.name2haplogroup.values());
	}

private HaploGroupImpl parse(HaploGroupImpl parent,Element root) {
	if(!root.hasAttribute("name"))  throw new IllegalArgumentException("expected @name");
	final String sn = root.getAttribute("name");
	if(StringUtils.isBlank(sn)) throw new IllegalArgumentException("expected non empty @name");
	if(this.name2haplogroup.containsKey(sn)) throw new IllegalArgumentException("duplicate name "+sn);
	final HaploGroupImpl hg = new HaploGroupImpl(sn);
	this.name2haplogroup.put(sn, hg);
	if(parent!=null) {
		parent.children.add(hg);
		hg.parent = parent;
		}
	for(Node c1=root.getFirstChild();c1!=null;c1=c1.getNextSibling()) {
		if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
		final Element E1= Element.class.cast(c1);
		if(E1.getNodeName().equals("haplogroup")) {
			this.parse(hg,E1);
			}
		else if(E1.getNodeName().equals("details")) {
			for(Node c2=E1.getFirstChild();c2!=null;c2=c2.getNextSibling()) {
				if(c2.getNodeType()!=Node.ELEMENT_NODE) continue;
				final Element E2= Element.class.cast(c2);
				if(E2.getNodeName().equals("poly")) {
					/* BUGGY not just position-base  e.g: '5899.XC'
					final String s=E2.getTextContent().trim();
					final int pos = Integer.parseInt(s.substring(0, s.length()-1));
					final Allele alt = Allele.create(s.substring( s.length()-1), false);
					hg.polymorphisms.add(new PolymorphismImpl(pos, alt));
					*/
					}
					
				}
			}
		}
	return hg;
	}

public boolean hasHaplogroup(final String sn) {
	return this.name2haplogroup.containsKey(sn);
	}

public HaploGroup getHaplogroupByName(final String sn) {
	return this.name2haplogroup.get(sn);
	}


public List<HaploGroup> getBaseHaplogroups() {
	return this.name2haplogroup.values().stream().
			filter(H->H.parent==null).
			collect(Collectors.toList());
	}


public static PhyloTree  load(final String urlOrFile) throws IOException {
	if(IOUtils.isRemoteURI(urlOrFile)) {
		return load(new InputSource(urlOrFile));
		}
	else
		{
		try(InputStream in= IOUtils.openPathForReading(Paths.get(urlOrFile))) {
			return load(new InputSource(in));
			}
		}
	}

public static PhyloTree  load(final InputSource inputSource) throws IOException {
	try {
		final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		dbf.setNamespaceAware(false);
		final DocumentBuilder db = dbf.newDocumentBuilder();
		final Document dom = db.parse(inputSource);
		final Element root= dom.getDocumentElement();
		if(root==null) throw new IOException("cannot find root element");
		if(!root.getNodeName().equals("phylotree")) throw new IOException("expected <phylotree> root but got "+root.getNodeName());
		final PhyloTree phyloTree = new PhyloTree();
		for(Node c1=root.getFirstChild();c1!=null;c1=c1.getNextSibling()) {
			if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
			Element E1= Element.class.cast(c1);
			if(!E1.getNodeName().equals("haplogroup")) continue;
			phyloTree.parse(null,E1);
			}
		return phyloTree;
		}
	catch(ParserConfigurationException|SAXException err) {
		throw new IOException(err);
		}
	}
}
