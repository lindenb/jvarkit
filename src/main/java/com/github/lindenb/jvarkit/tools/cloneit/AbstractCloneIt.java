package com.github.lindenb.jvarkit.tools.cloneit;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.github.lindenb.jvarkit.util.bio.Rebase;

import htsjdk.samtools.util.CoordMath;

public class AbstractCloneIt {

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



protected interface Enzyme
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

protected class EnzymeImpl implements Enzyme
	{
	final int _5_3;
	final Rebase.Enzyme rebaseEnz;
	EnzymeImpl(final Rebase.Enzyme rebaseEnz) {
		this.rebaseEnz = rebaseEnz;
		final String decl = rebaseEnz.getDecl();
		final int k = decl.indexOf('^');
		if( k == -1) throw new IllegalArgumentException(decl+" missing ^");
		this._5_3=k;//TODO
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


protected abstract class AbstractPlasmid 
	implements Plasmid
	{
	protected byte sequence[];
	protected String name;
	protected List<Site> sites = new ArrayList<>();
	protected Comparator<Site> my_comparator = (S1,S2) -> Integer.compare(S1.getPosition(), S2.getPosition());
	
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
		return false;
		}
	Plasmid digest(final List<Enzyme> rebase) {
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
		return this;
		}
	}

protected class VectorPlasmidImpl extends AbstractPlasmid implements VectorPlasmid
	{
	private Range polylinker;
	@Override
	public Range getPolylinker() {
		return this.polylinker;
		}
	}
protected class InsertPlasmidImpl extends AbstractPlasmid implements InsertPlasmid
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

}
