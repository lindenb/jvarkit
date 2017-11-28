package com.github.lindenb.jvarkit.tools.cloneit;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Algorithms;
import com.github.lindenb.jvarkit.util.bio.Rebase;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequence;
import com.github.lindenb.jvarkit.util.bio.fasta.FastaSequenceReader;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.sun.javafx.scene.traversal.Algorithm;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;

public abstract class AbstractCloneIt extends Launcher {

protected static class Range
	{
	final int start;
	final int end;
	Range(int start,int end) {
		this.start = start;
		this.end = end;
		}
	public int getStart() {
		return start;
		}
	public int getEnd() {
		return end;
	}
	public boolean overlaps(final Range other) {
		return CoordMath.overlaps(getStart(), getEnd(), other.getStart(),other.getEnd());
		}
	@Override
	public String toString() {
		return "["+this.start+"-"+this.end+"]";
		}
	}
	
protected interface Plasmid
	{
	public char charAt(int index);
	public int length();
	public String getName();
	public List<Site> getSites();
	public void digest(final List<Enzyme> rebase);
	}

protected interface InsertPlasmid extends Plasmid
	{
	public Range getPolylinker5();
	public Range getPolylinker3();
	}
protected interface VectorPlasmid extends Plasmid
	{
	public Range getPolylinker();
	
	public default List<Site> getSitesInPolylinker() {
		return this.getSites().stream().filter(
				S->true
				).collect(Collectors.toList());
		}
	}



protected static interface Enzyme
	{
	public String getName();
	public String getBases();
	public String getDeclaredBases();
	public int length();
	public char charAt(int index);
	public int get5_3();
	public default int get3_5() { return length()-get5_3();}
	public default int getOverhangLength() { return get5_3()-get3_5();}
	}

protected static class EnzymeImpl implements Enzyme
	{
	final int _5_3;
	final Rebase.Enzyme rebaseEnz;
	final int _hash;
	EnzymeImpl(final Rebase.Enzyme rebaseEnz) {
		this.rebaseEnz = rebaseEnz;
		final String decl = rebaseEnz.getDecl();
		final int k = decl.indexOf('^');
		if( k == -1) throw new IllegalArgumentException(decl+" missing ^");
		this._5_3=k;//TODO
		this._hash = decl.hashCode();
		}
	@Override
	public String getName() {
		return this.rebaseEnz.getName();
		}
	@Override
	public String getDeclaredBases() {
		return this.rebaseEnz.getDecl();
		}
	
	@Override
	public String getBases() {
		return this.rebaseEnz.getBases();
		}
	
	@Override
	public char charAt(int index) {
		return getBases().charAt(index);
		}
	@Override
	public int length() {
		return getBases().length();
		}
	@Override
	public int get5_3() {
		return _5_3;
		}
	@Override
	public int hashCode() {
		return this._hash;
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==null) return false;
		if(obj==this) return true;
		if(!(obj instanceof Enzyme)) return false;
		return this.getDeclaredBases().equals(Enzyme.class.cast(obj).getDeclaredBases());
		}
	@Override
	public String toString() {
		return rebaseEnz.toString();
		}
	}

protected interface Site
	{
	public Enzyme getEnzyme();
	public Plasmid getPlasmid();
	public int getPosition();
	public default int getMaxPosition() { return getPosition() + getEnzyme().length();}
	public default int get5_3() { return getPosition() + getEnzyme().get5_3();}
	public default int get3_5() { return getPosition() + getEnzyme().get3_5();}
	public default int getOverhangLength() { return get5_3()-get3_5();}
	public default boolean isStickyWith(final Site right) {
		if(this.getOverhangLength()!=right.getOverhangLength()) return false;
		for(int i= 0;i< this.getOverhangLength();i++)
			{
			char b1 = this.getPlasmid().charAt(this.get5_3()+i);
			char b2 = right.getPlasmid().charAt(right.get5_3()+i);
			if(Character.toUpperCase(b1)!=Character.toUpperCase(b2)) return false;
			}
		return true;
		}
	}

protected interface SiteWithTreatment extends Site
	{
	public boolean isTreated();
	public default int get5_3() { return isTreated()?0:getPosition() + getEnzyme().get5_3();}
	public default int get3_5() { return isTreated()?0:getPosition() + getEnzyme().get3_5();}
	}


protected static abstract class AbstractPlasmid 
	implements Plasmid
	{
	protected byte sequence[];
	protected String name;
	protected final List<Site> sites = new ArrayList<>();
	protected final Comparator<Site> my_comparator = (S1,S2) -> Integer.compare(S1.getPosition(), S2.getPosition());
	
	private static class CompositeList<S> extends AbstractList<S>
		{
		final List<S> L;
		final List<S> R;
		CompositeList(final List<S> L,final List<S> R) {
			this.L = L;
			this.R = R;
			}
		@Override
		public S get(int index) {
			return (index<L.size()?L.get(index):R.get(index-L.size()));
			}
		@Override
		public int size() {
			return this.L.size() + this.R.size();
			}
		}
	
	
	private class SiteImpl implements Site
		{
		final Enzyme enz;
		final int pos;
		SiteImpl(final Enzyme enz,int pos) {
			this.enz = enz;
			this.pos = pos;
			}
		@Override
		public Enzyme getEnzyme() {
			return this.enz;
			}
		@Override
		public int getPosition() {
			return this.pos;
			}
		@Override
		public Plasmid getPlasmid() {
			return AbstractPlasmid.this;
			}
		}

	
	@Override
	public String getName() {
		return name;
		}
	
	@Override
	public List<Site> getSites() {
		return sites;
		}
	
	@Override
	public int length() {
		return sequence.length;
		}
	
	@Override
	public char charAt(int index) {
		return (char)this.sequence[index];
		}
	private boolean match(char e, char p) {
		return Rebase.compatible(p, e);
		}
	@Override
	public void digest(final List<Enzyme> rebase) {
		this.sites.clear();
		for(int i=0;i< this.length();i++) {
			for(final Enzyme enz: rebase) {
				int x=0;
				while(x< enz.length()) {
					if(!match(enz.charAt(x),this.charAt(i+x))) break;
					++x;
					}
				if(x!=enz.length()) continue;
				this.sites.add(new SiteImpl(enz,i));
				}
			}
		
		Collections.sort(this.sites,this.my_comparator);
		}
	

	
	}

protected static class VectorPlasmidImpl extends AbstractPlasmid implements VectorPlasmid
	{
	private Range polylinker;
	@Override
	public Range getPolylinker() {
		return this.polylinker;
		}
	}
protected static class InsertPlasmidImpl extends AbstractPlasmid implements InsertPlasmid
	{
	private Range polylinker5;
	private Range polylinker3;
	@Override
	public Range getPolylinker5() {
		return this.polylinker5;
		}
	@Override
	public Range getPolylinker3() {
		return this.polylinker3;
		}
	}

protected class SiteWithTreatmentImpl implements SiteWithTreatment
	{
	final Site delegate;
	final boolean treated;
	SiteWithTreatmentImpl(final Site delegate, boolean treated) {
		this.delegate = delegate;
		this.treated = treated;
		}
	@Override
	public Enzyme getEnzyme() { return this.delegate.getEnzyme();}
	@Override
	public int getPosition() { return this.delegate.getPosition();}
	@Override
	public Plasmid getPlasmid() { return this.delegate.getPlasmid();}
	@Override
	public boolean isTreated() { return this.treated;}
	}

public static abstract class BasicPlasmidLoader
	{
	protected int[] parse_intervals(final String s,int length,int expect) {
		final String tokens[]=s.split("[,]");
		if(tokens.length!=expect)
			{
			throw new JvarkitException.UserError("expected "+expect+ " number in "+s);
			}
		int pos[]=new int[tokens.length];
		for(int i=0;i< tokens.length;++i) {
			final int v;
			try { v= Integer.parseInt(tokens[i]);}
			catch(NumberFormatException err) {
				throw new JvarkitException.UserError("bad coordinate in "+s);
				}
			if(v<1) throw new JvarkitException.UserError("bad coordinate "+v+"<1 in "+s);
			if(v>length) throw new JvarkitException.UserError("bad coordinate "+v+">"+length+" in "+s);
			if(i>0 && pos[i-1]>=v) throw new JvarkitException.UserError("bad coordinate "+v+" in "+s+" number should increase.");
			pos[i]=v-1;
		}
		return pos;
		}
	
	protected FastaSequence load(final File fastaFile) throws IOException
		{
		return new FastaSequenceReader().readOne(fastaFile);
		}
	}


public static class VectorPlasmidLoader
	extends BasicPlasmidLoader
	{
	@Parameter(names={"-vector","--vector"},description="Path to a fasta file containing one fasta file for the vector")
	private File fastaFile;
	@Parameter(names={"-vp","--vector-polylinker"},description="Vector polylinker. A group of 2 coordinates (1-based) for the begining and the inclusive end of the polylinker")
	private String polylinker="";
	
	VectorPlasmid create() throws IOException{
		FastaSequence f=this.load(this.fastaFile);
		int pos[]=this.parse_intervals(this.polylinker, f.length(), 2);
		VectorPlasmidImpl p=new VectorPlasmidImpl();
		p.polylinker=new Range(pos[0], pos[1]);
		p.sequence = f.toByteArray();
		p.name = f.getName();
		return p;
		}

	}

public static class InsertPlasmidLoader
	extends BasicPlasmidLoader
	{
	@Parameter(names={"-insert","--insert"},description="Path to a fasta file containing one fasta file for the insert")
	private File fastaFile;
	@Parameter(names={"-ip","--insert-polylinker"},description="Insert polylinker. A group of 4 coordinates (1-based) for the begining and the inclusive end of the polylinker")
	private String polylinker="";
	
	InsertPlasmid create() throws IOException{
		final FastaSequence f=this.load(this.fastaFile);
		final int pos[]=this.parse_intervals(this.polylinker, f.length(),4);
		final InsertPlasmidImpl p=new InsertPlasmidImpl();
		p.polylinker5=new Range(pos[0], pos[1]);
		p.polylinker3=new Range(pos[2], pos[3]);
		p.sequence = f.toByteArray();
		p.name = f.getName();
		return p;
		}
	}
private final Rebase rebase = Rebase.createDefaultRebase();

protected AbstractCloneIt() {
	}

protected Rebase getRebase() {
	return rebase;
	}


}
