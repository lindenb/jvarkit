/*
 * Disq
 *
 * MIT License
 *
 * Copyright (c) 2018-2019 Disq contributors
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.disq_bio.disq.impl.formats.bgzf;

import htsjdk.samtools.seekablestream.SeekableStream;
import java.io.Closeable;
import java.io.IOException;
import java.io.Serializable;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import com.github.lindenb.jvarkit.io.IOUtils;

public class BgzfBlockGuesser implements Closeable {

  protected static final int BGZF_MAGIC = 0x04088b1f;
  protected static final int BGZF_MAGIC_SUB = 0x00024342;
  protected static final int BGZF_SUB_SIZE = 4 + 2;

  protected SeekableStream in;
  protected String path;
  protected final ByteBuffer buf;

  public BgzfBlockGuesser(SeekableStream in, String path) {
    this.in = in;
    this.path = path;
    buf = ByteBuffer.allocate(8);
    buf.order(ByteOrder.LITTLE_ENDIAN);
  }

  public static class BgzfBlock implements Serializable {
    private static final long serialVersionUID = 1L;

    public String path;
    public long pos;
    public int cSize;
    public int uSize;
    public long end;
    private transient SeekableStream in;

    public BgzfBlock(String pa, long p, int cs, int us, long e, SeekableStream in) {
      path = pa;
      pos = p;
      cSize = cs;
      uSize = us;
      end = e;
      this.in = in;
    }

    /**
     * Signals that this is the last block in a partition, and that the guesser should close its
     * resources.
     *
     * @throws IOException if an I/O error occurs
     */
    public void end() throws IOException {
      in.close();
    }

    @Override
    public String toString() {
      return "BgzfBlock{"
          + "path="
          + path
          + ", pos="
          + pos
          + ", cSize="
          + cSize
          + ", uSize="
          + uSize
          + ", end="
          + end
          + '}';
    }
  }

  // Gives the compressed size on the side. Returns null if it doesn't find
  // anything.
  public BgzfBlock guessNextBGZFPos(long p, long end) {
    try {
      for (; ; ) {
        for (; ; ) {
          in.seek(p);
          IOUtils.readFully(in, buf.array(), 0, 4);
          int n = buf.getInt(0);

          if (n == BGZF_MAGIC) break;

          // Skip ahead a bit more than 1 byte if you can.
          if (n >>> 8 == BGZF_MAGIC << 8 >>> 8) ++p;
          else if (n >>> 16 == BGZF_MAGIC << 16 >>> 16) p += 2;
          else p += 3;

          if (p >= end) return null;
        }
        // Found what looks like a gzip block header: now get XLEN and
        // search for the BGZF subfield.
        final long p0 = p;
        p += 10;
        in.seek(p);
        IOUtils.readFully(in, buf.array(), 0, 2);
        p += 2;
        final int xlen = getUShort(0);
        final long subEnd = p + xlen;

        while (p < subEnd) {
          IOUtils.readFully(in, buf.array(), 0, 4);

          if (buf.getInt(0) != BGZF_MAGIC_SUB) {
            p += 4 + getUShort(2);
            in.seek(p);
            continue;
          }

          // Found it: this is close enough to a BGZF block, make it
          // our guess.

          // But find out the size before returning. First, grab bsize:
          // we'll need it later.
          IOUtils.readFully(in, buf.array(), 0, 2);
          int bsize = getUShort(0);

          // Then skip the rest of the subfields.
          p += BGZF_SUB_SIZE;
          while (p < subEnd) {
            in.seek(p);
            IOUtils.readFully(in, buf.array(), 0, 4);
            p += 4 + getUShort(2);
          }
          if (p != subEnd) {
            // Cancel our guess because the xlen field didn't match the
            // data.
            break;
          }

          // Now skip past the compressed data and the CRC-32.
          p += bsize - xlen - 19 + 4;
          in.seek(p);
          IOUtils.readFully(in, buf.array(), 0, 4);
          return new BgzfBlock(path, p0, (int) (p + 4 - p0), buf.getInt(0), end, in);
        }
        // No luck: look for the next gzip block header. Start right after
        // where we last saw the identifiers, although we could probably
        // safely skip further ahead. (If we find the correct one right
        // now, the previous block contained 0x1f8b0804 bytes of data: that
        // seems... unlikely.)
        p = p0 + 4;
      }
    } catch (IOException e) {
      return null;
    }
  }

  protected int getUShort(final int idx) {
    return (int) buf.getShort(idx) & 0xffff;
  }

  @Override
  public void close() {
    try {
      in.close();
    } catch (final IOException e) {
      throw new RuntimeException(e);
    }
  }
}
