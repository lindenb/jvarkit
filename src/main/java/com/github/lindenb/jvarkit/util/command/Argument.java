package com.github.lindenb.jvarkit.util.command;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Supplier;

public interface Argument<T> extends Supplier<T> {
	static final String HELP_SHORT_OPT="h";
	static final String HELP_LONG_OPT="help";
	static final String VERSION_LONG_OPT="version";

	
	public boolean isHidden();
	public boolean isDeprecated();
	public String getLongOpt();
	public String getShortOpt();
	public String getShortDesc();
	public String getLongDesc();
	public Map<String,String> getMetaData();
	public Set<String> getCategories();
	public boolean isNilleable();
	default public boolean isNull() { return get()==null;}
	public void set(final T v);
	default public boolean hasShortOpt() { return this.getShortOpt()!=null && !this.getShortOpt().isEmpty();}
	default public boolean hasLongOpt() { return this.getLongOpt()!=null && !this.getLongOpt().isEmpty();}
	public boolean canDecodeOption(final List<String> args,int optind);
	public void decodeOption(final List<String> args,int optind);
	public void postParsing();
	default public T getNonNull() {
		if(isNull()) throw new NullPointerException("Content is null");
		return this.get();
	}
}
