package com.github.lindenb.jvarkit.gatk;

import java.util.List;
import java.util.Optional;

public interface Gatk4Proxy {
	public void execute(final List<String> argv) throws Exception;
	public static Optional<Gatk4Proxy> getInstance() {
			try {
				Class<?> clazz=Class.forName("com.github.lindenb.jvarkit.gatk.Gatk4ProxyImpl");
				return Optional.of(Gatk4Proxy.class.cast(clazz.getDeclaredConstructor().newInstance()));
				}
			catch(Throwable err) {
				err.printStackTrace();
				return Optional.empty();
				}
		}
}
