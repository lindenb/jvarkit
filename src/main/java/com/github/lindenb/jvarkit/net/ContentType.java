/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.net;

import java.util.Optional;

/** classical content-types */
public class ContentType {
public static final ContentType TEXT_PLAIN = of("text/pain");
public static final ContentType TEXT_HTML = of("text/html");
public static final ContentType TEXT_XML = of("text/xml");
public static final ContentType TEXT_CSS = of("text/css");
public static final ContentType TEXT_JAVASCRIPT = of("text/javascript");
public static final ContentType IMAGE_SVG = of("image/svg+xml");
public static final ContentType APPLICATION_JSON = of("application/json");
public static final ContentType APPLICATION_OCTECT_STREAM = of("application/octet-stream");
public static final ContentType APPLICATION_ZIP = of("application/zip");
public static final ContentType APPLICATION_PDF = of("application/pdf");

public static Optional<ContentType> fromSuffix(String filename) {
	filename = filename.toLowerCase();
	if(filename.endsWith(".html")) return Optional.of(TEXT_HTML);
	if(filename.endsWith(".txt")) return Optional.of(TEXT_PLAIN);
	if(filename.endsWith(".json")) return Optional.of(APPLICATION_JSON);
	if(filename.endsWith(".xml")) return Optional.of(TEXT_XML);
	if(filename.endsWith(".js")) return Optional.of(TEXT_JAVASCRIPT);
	if(filename.endsWith(".css")) return Optional.of(TEXT_CSS);
	if(filename.endsWith(".svg")) return Optional.of(IMAGE_SVG);
	if(filename.endsWith(".zip")) return Optional.of(APPLICATION_ZIP);
	if(filename.endsWith(".pdf")) return Optional.of(APPLICATION_PDF);
	return Optional.empty();
	}
private final String mimeType;
private ContentType(final String mimeType) {
	this.mimeType=mimeType;
	}
private static ContentType of(final String mimeType) {
	return new ContentType(mimeType);
	}
public String getMimeType() {
	return mimeType;
	}
@Override
public String toString() {
	return getMimeType();
	}
}
