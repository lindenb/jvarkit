/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util.bio;

import java.util.Set;

import htsjdk.variant.variantcontext.Allele;

public class AcidNucleics {

/** return true if a is not null, not symbolic, not span not empty and all bases match 'isATGC */
public static boolean isATGC(final Allele a) {
	if(a==null || a.isBreakpoint() || a.isNoCall() || a.isSymbolic() || a.isSingleBreakend() || a.length()==0) return false;
	return isATGC(a.getDisplayString());
	}
	
/** return true if s is not null, not empty and all bases match 'isATGC */
public static boolean isATGC(final CharSequence c) {
	if(c==null || c.length()==0) return false;
	for(int i=0;i< c.length();i++) if(!isATGC(c.charAt(i))) return false;
	return true;
	}
	
/** return true is base is [A,T,G,C] regardless of the case */
public static boolean isATGC(char c) {
	switch(c) {
		case 'A': case 'a':
		case 'C': case 'c':
		case 'G': case 'g':
		case 'T': case 't':
			return true;
		default: return false;
		}
	}
public static boolean isATGC(byte c) {
	return isATGC((char)c);
	}

/** returns true of base in AGR */
public static boolean isPuryn(char iupac) {
	switch(iupac) {
	case 'A': case 'a' :
	case 'G': case 'g':
	case 'R': case 'r':
		return true;
	default: return false;
	}
}

/** returns true of base in TUCY */
public static boolean isPyrimidic(char iupac) {
	switch(iupac) {
	case 'T': case 't' :
	case 'U': case 'u':
	case 'C': case 'c':
	case 'Y': case 'y':
		return true;
	default: return false;
	}
}

/** return true if s is not null, not empty and all bases match 'isATGCN */
public static boolean isATGCN(final CharSequence c) {
	if(c==null || c.length()==0) return false;
	for(int i=0;i< c.length();i++) if(!isATGCN(c.charAt(i))) return false;
	return true;
	}
	
/** return true is base is [A,T,G,C,N] regardless of the case */
public static boolean isATGCN(char c) {
	switch(c) {
		case 'n': case 'N': return true;
		default: return isATGC(c);
		}
	}
public static boolean isATGCN(byte c) {
	return isATGCN((char)c);
	}

/** return true if a is not null, not symbolic, not span not empty and all bases match 'isATGCN */
public static boolean isATGCN(final Allele a) {
	if(a==null || a.isBreakpoint() || a.isNoCall() || a.isSymbolic() || a.isSingleBreakend() || a.length()==0) return false;
	return isATGCN(a.getDisplayString());
	}


/** return the reverse complement of the sequence */	
public static String reverseComplement(final CharSequence seq)
	{
	final StringBuilder b=new StringBuilder(seq.length());
	for(int i=0;i< seq.length();++i)
		{
		b.append(complement(seq.charAt((seq.length()-1)-i)));
		}
	return b.toString();
	}

/** returns the reverse complement of the character */	
public static char complement(char c)
	{
	switch(c)
		{
        case 'A': return 'T';
        case 'T': case 'U': return 'A';
        case 'G': return 'C';
        case 'C': return 'G';

        case 'a': return 't';
        case 't': case 'u': return 'a';
        case 'g': return 'c';
        case 'c': return 'g';

        case 'w': return 'w';
        case 'W': return 'W';

        case 's': return 's';
        case 'S': return 'S';

        case 'y': return 'r';
        case 'Y': return 'R';

        case 'r': return 'y';
        case 'R': return 'Y';

        case 'k': return 'm';
        case 'K': return 'M';

        case 'm': return 'k';
        case 'M': return 'K';

        case 'b': return 'v';
        case 'd': return 'h';
        case 'h': return 'd';
        case 'v': return 'b';


        case 'B': return 'V';
        case 'D': return 'H';
        case 'H': return 'D';
        case 'V': return 'B';

        case 'N': return 'N';
        case 'n': return 'n';
		default: return 'N';
		}
	}

/** return the 'weight' of a degenerate base e.g: A=1, N=0, W=0.5, B=0.25, used with restriction enzymes */
public static float weight(char c)
	{
	switch(c)
		{
        case 'A': 
        case 'T':
        case 'U': 
        case 'G': 
        case 'C': 
        case 'a':
        case 't':
        case 'u': 
        case 'g': 
        case 'c': return 1f;

        case 'w': 
        case 'W': 

        case 's': 
        case 'S': 

        case 'y': 
        case 'Y': 

        case 'r': 
        case 'R':

        case 'k':
        case 'K': 

        case 'm': 
        case 'M': return 0.5f;

        case 'b': 
        case 'd': 
        case 'h': 
        case 'v': 
        case 'B': 
        case 'D': 
        case 'H': 
        case 'V': return 0.25f;

        case 'N':
        case 'n':
		default: return 0.f;
		}
	}
/** convert degenerate base to array of ATGC */
public static char[] degenerateToBases(final char base) {
	switch(base) {
	    case 'A':
	    case 'C':
	    case 'G':
	    case 'T':
	    case 'U': return new char[]{base};
	    case 'R': return new char[]{'A','G'}; 
	    case 'Y': return new char[]{'C','T'}; 
	    case 'S': return new char[]{'G','C'};
	    case 'W': return new char[]{'A','T'};
	    case 'K': return new char[]{'G','T'};
	    case 'M' :return new char[]{'A','C'};
	    case 'B': return new char[]{'C','G','T'};
	    case 'D': return new char[]{'A','G','T'};
	    case 'H': return new char[]{'A','C','T'};
	    case 'V': return new char[]{'A','C','G'};
	    case 'N': return new char[]{'A','C','G','T'};
	    //
	    case 'a':
	    case 'c':
	    case 'g':
	    case 't':
	    case 'u': return new char[]{base};
	    case 'r': return new char[]{'a','g'}; 
	    case 'y': return new char[]{'c','t'}; 
	    case 's': return new char[]{'g','c'};
	    case 'w': return new char[]{'a','t'};
	    case 'k': return new char[]{'g','t'};
	    case 'm' :return new char[]{'a','c'};
	    case 'b': return new char[]{'c','g','t'};
	    case 'd': return new char[]{'a','g','t'};
	    case 'h': return new char[]{'a','c','t'};
	    case 'v': return new char[]{'a','c','g'};
	    case 'n': return new char[]{'a','c','g','t'};
	    default: throw new IllegalArgumentException("bad DNA base:"+base);
		}
	}

/** looks like a nucleotide in the IUPAC classification */
public static boolean isIUPAC(final char base) {
	switch(base) {
	    case 'A':
	    case 'C':
	    case 'G':
	    case 'T':
	    case 'U':
	    case 'R':
	    case 'Y':
	    case 'S':
	    case 'W':
	    case 'K':
	    case 'M':
	    case 'B':
	    case 'D':
	    case 'H':
	    case 'V':
	    case 'N':
	    //
	    case 'a':
	    case 'c':
	    case 'g':
	    case 't':
	    case 'u':
	    case 'r':
	    case 'y':
	    case 's':
	    case 'w':
	    case 'k':
	    case 'm':
	    case 'b':
	    case 'd':
	    case 'h':
	    case 'v':
	    case 'n': return true;
	    default: return false;
		}
	}


/** TODO: not tested */
public static char basesToDegenerate(final Set<Character> bases) {
	switch(bases.size())
		{
		case 0 : throw new IllegalArgumentException("empty set of bases");
		case 1 :
			{
			final char c = bases.iterator().next();
			switch(c) {
				case 'u': case 'U':
				case 'a': case 'c': case 'g': case 't':
				case 'A': case 'C': case 'G': case 'T': return c;
				default : throw new IllegalArgumentException("Not ATGC: "+c);
				}
			}
		case 2:
			{
			if( (bases.contains('a') || bases.contains('A')) &&
			    (bases.contains('g') || bases.contains('G'))
				) return 'R';
			if( (bases.contains('c') || bases.contains('C')) &&
				(bases.contains('t') || bases.contains('T') || bases.contains('u') || bases.contains('U'))
				) return 'Y';
			if( (bases.contains('c') || bases.contains('C')) &&
				(bases.contains('g') || bases.contains('G'))
				) return 'S';
			if( (bases.contains('a') || bases.contains('A')) &&
				(bases.contains('t') || bases.contains('T') || bases.contains('u') || bases.contains('U'))
				) return 'W';
			if( (bases.contains('g') || bases.contains('G')) &&
				(bases.contains('t') || bases.contains('T') || bases.contains('u') || bases.contains('U'))
				) return 'K';
			if( (bases.contains('a') || bases.contains('A')) &&
				(bases.contains('c') || bases.contains('C'))
				) return 'M';
			throw new IllegalArgumentException("bad couple of bases "+bases);
			}
		case 3:
			{
			if( (bases.contains('c') || bases.contains('C')) &&
				(bases.contains('g') || bases.contains('G')) &&
				(bases.contains('t') || bases.contains('T') || bases.contains('u') || bases.contains('U'))
				) return 'B';
			if( (bases.contains('a') || bases.contains('A')) &&
				(bases.contains('g') || bases.contains('G')) &&
				(bases.contains('t') || bases.contains('T') || bases.contains('u') || bases.contains('U'))
				) return 'D';
			if( (bases.contains('a') || bases.contains('A')) &&
				(bases.contains('c') || bases.contains('C')) &&
				(bases.contains('t') || bases.contains('T') || bases.contains('u') || bases.contains('U'))
				) return 'H';
			if( (bases.contains('a') || bases.contains('A')) &&
				(bases.contains('c') || bases.contains('C')) &&
				(bases.contains('g') || bases.contains('G'))
				) return 'V';
			throw new IllegalArgumentException("bad triple of bases "+bases);
			}
		case 4:
			{
			if( (bases.contains('a') || bases.contains('A')) &&
				(bases.contains('c') || bases.contains('C')) &&
				(bases.contains('g') || bases.contains('G')) &&
				(bases.contains('t') || bases.contains('T') || bases.contains('u') || bases.contains('U'))
				) return 'N';
			throw new IllegalArgumentException("bad group of 4 bases "+bases);
			}
		default : throw new IllegalArgumentException("bad number of bases");
		}
	}

	
}
