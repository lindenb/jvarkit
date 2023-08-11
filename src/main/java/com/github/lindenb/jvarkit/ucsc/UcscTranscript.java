/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.ucsc;

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.bed.BedInterval;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript.CodingRNA;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.KozakSequence;

import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;

/**
 * Describe a UCSC transcript in the genepred format
 * @author lindenb
 *
 */
public interface UcscTranscript extends Feature , BedInterval{
/** get Transcript Id */
public String getTranscriptId();
/** get name2 or 'gene' if available */
public String getName2();
public Strand getStrand();
public boolean isPositiveStrand();
public boolean isNegativeStrand();
public int getTxStart();
public int getTxEnd();

/** return CDS start or throw an exception if it is not a protein coding transcript */
public int getCdsStart();
/** return CDS end or throw an exception if it is not a protein coding transcript */
public int getCdsEnd();


/** return sum of length on reference of exons */
public default int getTranscriptLength() {
	return getExons().stream().mapToInt(E->E.getLengthOnReference()).sum();
	}

/** get TxStart */
@Override
public default int getStart() {
	return getBedStart() +1 ;
	}

/** get TxEnd */
@Override
public default  int getEnd() {
	return getBedEnd();
	}
/** get count of exons */
public int getExonCount();

public default boolean hasIntrons() {
	return this.getExonCount() > 1;
}

public default int getIntronCount() {
	return Math.max(0, getExonCount()-1);
	}

public int getExonStart(int idx);
public int getExonEnd(int idx);


public boolean isProteinCoding();

public default Intron getIntron(int idx) {
	return new Intron(this,idx);
	}

public default Exon getExon(int idx) {
	return new Exon(this,idx);
	}

public default List<Exon> getExons() {
	return new AbstractList<Exon>() {
		@Override
		public Exon get(int idx) {
			return getExon(idx);
			}
		@Override
		public int size() {
			return getExonCount();
			}
		};
	}



public default List<Intron> getIntrons() {
	if(!hasIntrons()) return Collections.emptyList();
	return new AbstractList<Intron>() {
		@Override
		public Intron get(int idx) {
			return getIntron(idx);
			}
		@Override
		public int size() {
			return getIntronCount();
			}
		};
	}

public default List<CDS> getCDSs() {
	if(!isProteinCoding()) Collections.emptyList();
	return getExons().stream().
		filter(EX->CoordMath.overlaps(
				EX.getStart(), EX.getEnd(),
				getCdsStart(), getCdsEnd()
				)).
		map(EX->new CDS(EX)).
		collect(Collectors.toList());
	}




public abstract class Component implements Locatable, BedInterval {
	public abstract UcscTranscript getTranscript();
	
	@Override
	public int hashCode() {
		int i = getContig().hashCode();
		i= i*31 + getTranscript().getTranscriptId().hashCode();
		i= i*31 + Integer.hashCode(getStart());
		i= i*31 + Integer.hashCode(getEnd());
		return i;
		}
	
	
	public boolean containsGenomicLoc1(int pos1) {
		return getStart() <= pos1 && pos1 <= getEnd();
	}
	
	public abstract String getName();
	
	@Override
	public final int getBedStart() {
		return getStart()-1;
		}
	
	@Override
	public final int getBedEnd() {
		return getEnd();
		}
	
	public Strand getStrand() {
		return getTranscript().getStrand();
		}
	
	public boolean isPositiveStrand() {
		return getTranscript().isPositiveStrand();
		}
	
	public boolean isNegativeStrand() {
		return getTranscript().isNegativeStrand();
		}
	
	@Override
	public String getContig() {
		return getTranscript().getContig();
		}
	@Override
	public String toString() {
		return getClass().getSimpleName()+":"+getTranscript().getTranscriptId()+":"+getContig()+":"+getStart()+"-"+getEnd();
		}
	
	public Interval toInterval() {
		return new Interval(getContig(), getStart(), getEnd(), isNegativeStrand(), getName());
		}
	}



/************************************************************************
 * 
 * Exon
 *
 */
public class Exon extends Component {
	private final UcscTranscript owner;
	private final int exon_index;
	Exon(final UcscTranscript owner,int exon_index) {
		this.owner = owner;
		this.exon_index = exon_index;
		}
	
	@Override
	public UcscTranscript getTranscript() {
		return this.owner;
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(this==obj) return true;
		if(obj==null || !(obj instanceof Exon)) return false;
		final Exon o = Exon.class.cast(obj);
		return exon_index == o.exon_index &&
				contigsMatch(o) &&
				getTranscript().getTranscriptId().equals(o.getTranscript().getTranscriptId());
		}
	
	/** return 1-based exon index , according to strand */
	public int getHumanIndex() {
		return isPositiveStrand()?exon_index+1:getTranscript().getExonCount()-exon_index;
		}
	
	
	
	@Override
	public int getStart() {
		return getTranscript().getExonStart(this.exon_index);
		}
	@Override
	public int getEnd() {
		return getTranscript().getExonEnd(this.exon_index);
		}
	@Override
	public String getName() {
		return getTranscript().getTranscriptId()+".Exon"+ getHumanIndex();
		}
	}



/************************************************************************
 * 
 * Intron
 *
 */
public class Intron extends Component {
	private final UcscTranscript owner;
	private final int intron_index;
	private Intron(final UcscTranscript owner,int intron_index) {
		this.owner = owner;
		this.intron_index = intron_index;
		}
	
	@Override
	public UcscTranscript getTranscript() {
		return this.owner;
		}

	/** return 1-based intron index , according to strand */
	public int getHumanIndex() {
		return isPositiveStrand()?intron_index+1:getTranscript().getIntronCount()-intron_index;
		}
	
	@Override
	public int getStart() {
		return getTranscript().getExonEnd(this.intron_index) +1;
		}
	@Override
	public int getEnd() {
		return getTranscript().getExonStart(this.intron_index+1) -1;
		}
	@Override
	public String getName() {
		return getTranscript().getTranscriptId()+".Intron"+getHumanIndex();
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(this==obj) return true;
		if(obj==null || !(obj instanceof Intron)) return false;
		final Intron o = Intron.class.cast(obj);
		return intron_index==o.intron_index && 
				contigsMatch(o) && 
				getTranscript().getTranscriptId().equals(o.getTranscript().getTranscriptId());
		}
	/** get Surrounding exons */
	public Exon[] getEnclosingExons() {
		return new Exon[] {
				getTranscript().getExon(this.intron_index),
				getTranscript().getExon(this.intron_index+1),
			};
		}
	}

/** ===================================================================*/
public abstract class ExonComponent extends Component {
	private final Exon exon;
	private ExonComponent(final Exon exon) {
		this.exon = exon;
		}
	
	@Override
	public UcscTranscript getTranscript() {
		return this.exon.getTranscript();
		}
	
	public final Exon getExon() {
		return this.exon;
		}
	@Override
	public String getName() {
		return getTranscript().getTranscriptId()+"."+getClass().getSimpleName() + "."+getStart()+"-"+getEnd();
		}
	}



public class CDS extends ExonComponent {
	private CDS(final Exon exon) {
		super(exon);
		}
	@Override
	public int getStart() {
		return Math.max(
			getExon().getStart(),
			getTranscript().getCdsStart()
			);
		}
	@Override
	public int getEnd() {
		return Math.min(
			getExon().getEnd(),
			getTranscript().getCdsEnd()
			);
		}
	@Override
	public boolean equals(Object obj) {
		if(this==obj) return true;
		if(obj==null || !(obj instanceof CDS)) return false;
		final CDS o = CDS.class.cast(obj);
		return getStart()==o.getStart() && 
				getEnd()==o.getEnd() && 
				contigsMatch(o) && 
				getTranscript().getTranscriptId().equals(o.getTranscript().getTranscriptId());
		}
	}


public default boolean hasUTR5() {
	if(!isProteinCoding()) return false;
	if(isPositiveStrand()) {
		return getTxStart() < getCdsStart();
		}
	else
		{
		return getCdsEnd() < getTxEnd();
		}
	}



public default List<UTR5> getUTR5() {
	if(!hasUTR5()) return Collections.emptyList();
	final List<UTR5> L = new ArrayList<>();
	if(isPositiveStrand()) {
		for(int i=0;i< getExonCount();i++) {
			final Exon ex = getExon(i);
			if(ex.getStart() >= getCdsStart()) break;
			L.add(new UTR5(ex));
			}
		}
	else
		{
		for(int i=0;i< getExonCount();i++) {
			final Exon ex = getExon(i);
			if(ex.getEnd() <= getCdsEnd()) continue;
			L.add(new UTR5(ex));
			}
		}
	return L;
	}

/** ===================================================================*/
public abstract class UTR extends ExonComponent {
	protected UTR(final Exon exon) {
		super(exon);
		}
	}
/** ===================================================================*/
public class UTR5 extends UTR {
	private UTR5(final Exon exon) {
		super(exon);
		}
	@Override
	public int getStart() {
		if(isPositiveStrand()) {
			return getExon().getStart();
			}
		else
			{
			return Math.max(
					getExon().getEnd(),
					getTranscript().getCdsEnd()
					);
			}
		}
	@Override
	public int getEnd() {
		if(isPositiveStrand()) {
			return Math.min(
				getExon().getEnd(),
				getTranscript().getCdsStart() -1
				);
			}
		else
			{
			return getExon().getEnd();
			}
		}
	}
/** ===================================================================*/
public class UTR3 extends UTR {
	private UTR3(final Exon exon) {
		super(exon);
		}
	@Override
	public int getStart() {
		if(isPositiveStrand()) {
			return Math.max(
				getExon().getStart(),
				getTranscript().getCdsEnd()
				);
			}
		else
			{
			return getExon().getStart();
			}
		}
	@Override
	public int getEnd() {
		if(isPositiveStrand()) {
			return getExon().getEnd();
			}
		else
			{
			return Math.min(
				getExon().getEnd(),
				getTranscript().getCdsStart() -1
				);
			}
		}

	}
/** ===================================================================*/
/** ===================================================================*/
public abstract class Codon extends Component {
	private final UcscTranscript owner;
	private final int indexes[]=new int[3];
	protected Codon(final UcscTranscript owner,int pos1) {
		this.owner = owner;
		if(!owner.isProteinCoding()) throw new IllegalArgumentException("Not a protein coding transcript");
		if(isPositiveStrand()) {
			for(int i=0;i< owner.getExonCount();i++) {
				Exon ex = owner.getExon(i);
				if(!ex.containsGenomicLoc1(pos1)) continue;
				indexes[0] = pos1;
				
				}
			}
		else
			{
			
			}
		}
	@Override
	public UcscTranscript getTranscript() {
		return owner;
		}
	@Override
	public int getStart() {
		return 0;
		}
	@Override
	public int getEnd() {
		return 0;
		}
	}
/** ===================================================================*/
public class StartCodon extends Codon {
	private StartCodon(final UcscTranscript owner) {
		super(owner,owner.getCdsStart());
		}
	@Override
	public String getName() {
		return getTranscript().getTranscriptId()+".start_codon";
		}
	}
/** ===================================================================*/
public class StopCodon extends Codon {
	private StopCodon(final UcscTranscript owner) {
		super(owner,owner.getCdsEnd());
		}
	
	@Override
	public String getName() {
		return getTranscript().getTranscriptId()+".stop_codon";
		}
	}

public interface RNA extends CharSequence, Locatable {
	public UcscTranscript getTranscript();
	/** get a description for this RNA */
	public String getDescription();
	public default String getContig() { return getTranscript().getContig(); }
	public default boolean isPositiveStrand() { return getTranscript().isPositiveStrand();}
	public default boolean isNegativeStrand() { return getTranscript().isNegativeStrand();}
	
	/** return true if codon starting a pos0 is an ATG */
	public boolean isATG(int pos0);
	/** return true if codon starting a pos0 is a STOP */
	public boolean isStop(int pos0);
}

public interface MessengerRNA extends RNA {
	public default boolean isProteinCoding() {
		return getTranscript().isProteinCoding();
		}
	public CodingRNA getCodingRNA();
}


public interface CodingRNA extends RNA {
	public Peptide getPeptide();
	public KozakSequence getKozakSequence();
	public MessengerRNA getMessengerRNA();
}


public interface FivePrimeUpstreamORF extends RNA {
}


/** return mRNA for this transcript */
public MessengerRNA getMessengerRNA(CharSequence chromosomeSequence);


/** return coding RNA for this transcript */
public default CodingRNA getCodingRNA(CharSequence chromosomeSequence) {
	return getMessengerRNA(chromosomeSequence).getCodingRNA();
	}

public interface Peptide extends CharSequence {
	public CodingRNA getCodingRNA();
	public int[] convertToGenomicCoordinates(int pepPos0);
}

}
