/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.table;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

public abstract class AbstractTableExporter {
private boolean printEmptyTable = false;
private boolean removeEmptyColumns = true;

public void setPrintEmptyTable(boolean printEmptyTable) {
	this.printEmptyTable = printEmptyTable;
}
public boolean isPrintEmptyTable() {
	return printEmptyTable;
	}
public void setRemoveEmptyColumns(boolean removeEmptyColumns) {
	this.removeEmptyColumns = removeEmptyColumns;
}
public boolean isRemoveEmptyColumns() {
	return removeEmptyColumns;
	}
public abstract void saveTableTo(final Table table,final PrintWriter pw) throws IOException;

public void saveTableTo(final Table table,final File file) throws IOException {
	try(final PrintWriter pw=new PrintWriter(file)) {
		this.saveTableTo(table,pw);
		pw.flush();
		pw.close();
		}
	}
}
