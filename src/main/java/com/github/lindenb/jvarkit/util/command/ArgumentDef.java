package com.github.lindenb.jvarkit.util.command;

import java.io.PrintWriter;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.ParameterizedType;
import java.lang.reflect.Type;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Function;
import java.util.regex.Pattern;

import org.apache.commons.cli.Option;

import htsjdk.samtools.util.StringUtil;


public class ArgumentDef
	{
	private final Object object;
	private final Field field;
	private final Argument argument;
	private final Option option ;
	public boolean visited=false;
	ArgumentDef(final Object object,final Field field,final Argument argument) {
		this.object = object;
		this.field = field;
		this.argument = argument;
		if(this.argument.opt()==null || this.argument.opt().isEmpty()) {
			throw new IllegalArgumentException("in "+ field +" opt cannot be empty");
			}
		if(this.argument.opt().contains("-")) {
			throw new IllegalArgumentException("in "+ field +" opt cannot contains hyphen "+this.argument.opt());
			}
		if(this.argument.opt().equals("h")) {
			throw new IllegalArgumentException("in "+ field +" opt -h is reserved");
			}
		final Option.Builder builder = Option.builder(this.getOpt());
		if(!StringUtil.isBlank(this.argument.longopt())) builder.longOpt(this.argument.longopt());
		builder.desc(this.argument.description());
		if( this.hasArgument()) {
			builder.argName(StringUtil.isBlank(this.argument.argname())?this.argument.argname():this.getOpt().toUpperCase());
			
			if(this.isMultiple()) {
				builder.hasArgs();
				}
			else
				{
				builder.hasArg(true);
				}
			}
		else
			{
			builder.hasArg(false);
			}
		
		this.option = builder.build();
		}
	
	public Field getField() {
		return this.field;
		}
	
	public Option getOption() {
		return this.option;
		}
	
	public String getOptDash() {
		return "-"+getOpt();
		}
	
	public String getOpt() {
		return argument.opt();
		}
	
	public boolean hasLongOpt() {
		if(!StringUtil.isBlank(this.argument.longopt()))  return true;
		return !argument.longopt().isEmpty();
		}

	public String getLongOptDash() {
		return getLongOpt()==null?null:"--"+getLongOpt();
		}
	
	public String getLongOpt() {
		if(!StringUtil.isBlank(this.argument.longopt())) return argument.longopt();
		return null;
		}
	
	void printHelp(PrintWriter pw) {
		if(this.argument.hidden()) return;
		StringBuilder sb=new StringBuilder(" ");
		sb.append(getOptDash());
			
		if(hasLongOpt()) {
			sb.append("|");
			sb.append(getLongOptDash());
			}
		sb.append(" ");
		if(hasArgument()) {
			sb.append("<").append(argument.argname()).append(">");
			}
		while(sb.length()%10!=0) sb.append(" ");
		sb.append("description");
		if(this.isMultiple()) {
			sb.append(" this argument can be invoked multiple times.");
			}
		else
			{
			try {
				sb.append(" Default: \"" +this.field.get(this.object)).append("\".");
				}
			catch(Exception err) {
			
				}
			}
		}
	
	public boolean isBoolean()
		{
		return Boolean.class==getField().getType() || Boolean.TYPE==getField().getType() ;
		}
	
	public boolean isMultiple() {
		final Type t=field.getGenericType();
		if(!( t instanceof ParameterizedType)) {
			return false;
			}
		ParameterizedType pt = ParameterizedType.class.cast(t);
		if(pt.getRawType() instanceof Class )
			{
			Class<?> clazz= (Class<?>)pt.getRawType();
			if(pt.getActualTypeArguments().length==1 &&  clazz.isAssignableFrom(Collection.class)) {
				return true;
				}
			}
		return false;
		}
	
	private Class<?> getDataClass()
		{
		final Type t=field.getGenericType();
		if(!( t instanceof ParameterizedType)) {
			return (Class<?>)t;
			}
		final ParameterizedType pt = ParameterizedType.class.cast(t);			
		return (Class<?>) pt.getActualTypeArguments()[0];	
		}
	
	private Function<String,Object> getStringConverter() {
		if(!hasArgument()) throw new IllegalStateException("Boolean can't have an argument");
		final Class<?> clazz = getDataClass();
		try {
			final Constructor<?> ctor =  clazz.getConstructor(String.class);
			return new Function<String, Object>()
				{
				@Override
				public Object apply(String t)
					{
					try {
						return ctor.newInstance(t);
						}
					catch(Exception err) {
						throw new RuntimeException("Cannot create a new instance of "+clazz+" from "+t);
						}
					}
				};
			}
		catch(final Exception err)
			{
			throw new RuntimeException("Cannot find a String constructor for class "+clazz,err);
			}
		}
	
	
	/** shortcut to !isBoolean */
	public boolean hasArgument() {
		return !isBoolean();
		}
	
	public boolean isHidden() {
		return this.argument.hidden();
		}
	
	Object parseValue(final String input) {
		if(!StringUtil.isBlank(this.argument.regex())) {
			final Pattern pat = Pattern.compile(this.argument.regex());
			if(!pat.matcher(input).matches()) {
				throw new IllegalArgumentException("Option "+getOptDash()+" should match regular expression : "+pat);
				}
			}
		try
			{
			return getStringConverter().apply(input);
			}
		catch (Exception e)
			{
			throw new RuntimeException(e);
			}
		}
	
	@Override
	public String toString()
		{
		StringBuilder sb=new  StringBuilder();
		sb.append(" has.argument=").append(hasArgument());
		sb.append(" is.boolean=").append(isBoolean());
		sb.append(" is.multiple=").append(isMultiple());
		sb.append(" opt=").append(getOptDash());
		if(hasLongOpt()) sb.append(" longopt=").append(getLongOptDash());
		sb.append(" class=").append(this.field.getDeclaringClass());
		sb.append(" field=").append(this.field.getName());
		sb.append(" type=").append(this.field.getType());
		sb.append(" type.primitive=").append(this.field.getType().isPrimitive());
		sb.append(" type=").append(this.field.getGenericType());
		return sb.toString();
		}
	
	}
