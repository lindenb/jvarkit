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
}
