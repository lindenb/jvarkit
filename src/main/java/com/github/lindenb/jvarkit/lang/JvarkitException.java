package com.github.lindenb.jvarkit.lang;

public class JvarkitException   {
@SuppressWarnings("serial")
public static class CommandLineError extends Error
	{
	CommandLineError(final String msg) {
		super(msg);
		}
	}
}
