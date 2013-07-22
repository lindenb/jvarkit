package com.github.lindenb.jvarkit.util.bio;

public class AcidNucleics {
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
		case 'N': return c;
		case 'n': return c;
		default: return 'N';
		}
	}
}
