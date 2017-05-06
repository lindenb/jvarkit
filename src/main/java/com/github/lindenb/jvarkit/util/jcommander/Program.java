package com.github.lindenb.jvarkit.util.jcommander;

import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.Target;

@Retention(java.lang.annotation.RetentionPolicy.RUNTIME)
@Target({ ElementType.TYPE })
public @interface Program {
	public String name() default "";
	public String description() default "";
	public String deprecatedMsg() default "";
	public String[] keywords() default {};
	public String[] authors() default "Pierre Lindenbaum PhD @yokofakun";
	public int[] biostars() default {};
	/** shall we generate the Markdown Documentation ? */
	public boolean generate_doc() default true;
	public com.github.lindenb.semontology.Term[] terms() default {};
}
