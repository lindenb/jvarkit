package com.github.lindenb.jvarkit.util.picard;

public class PicardException extends RuntimeException
{
	private static final long serialVersionUID = 1L;

	public PicardException(final String message) {
        super(message);
    }

    public PicardException(final String message, final Throwable throwable) {
        super(message, throwable);
    }

}
