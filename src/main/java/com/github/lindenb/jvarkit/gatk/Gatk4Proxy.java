package com.github.lindenb.jvarkit.gatk;

import java.util.List;
import java.util.Optional;

public interface Gatk4Proxy {
	public static final String MSG_COMPILATION="gatk4 engine is not available."
			+ " Not all versions of jvarkit will work ."
			+ "Jvarkit must be compiled with the java code of gatk4 using a similar compilation command: "
			+ "./gradlew jvarkit -Dgatk4.local.jar=/path/to/gatk/gatk-4.*/gatk-package-4.*-local.jar\n"
			+ "";
	public void execute(final List<String> argv) throws Exception;
	public static Optional<Gatk4Proxy> getInstance() {
			try {
				Class<?> clazz=Class.forName("com.github.lindenb.jvarkit.gatk.Gatk4ProxyImpl");
				return Optional.of(Gatk4Proxy.class.cast(clazz.getDeclaredConstructor().newInstance()));
				}
			catch(Throwable err) {
				return Optional.empty();
				}
		}
}
