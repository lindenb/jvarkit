package com.github.lindenb.jvarkit.util.bio;

import com.github.lindenb.jvarkit.lang.DelegateCharSequence;

public class RevCompCharSequence extends DelegateCharSequence
	{
	public RevCompCharSequence(CharSequence delegate) {
		super(delegate);
		}
	@Override
	public char charAt(int index)
		{
		return AcidNucleics.complement( super.charAt((super.length()-1)-index));
		}
	}
