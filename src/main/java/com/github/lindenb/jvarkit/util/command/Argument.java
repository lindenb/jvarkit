package com.github.lindenb.jvarkit.util.command;
import static java.lang.annotation.ElementType.FIELD;
import static java.lang.annotation.RetentionPolicy.RUNTIME;

import java.lang.annotation.Documented;
import java.lang.annotation.Inherited;
import java.lang.annotation.Retention;
import java.lang.annotation.Target;

@Retention(RUNTIME)
@Target(FIELD)
@Inherited
@Documented
public @interface Argument {
	public String opt();
	public String longopt() default "";
	public String description() default "";
	public String regex() default "";
	public String category() default "";
	public boolean hidden() default false;
}
