/*
 * Copyright 2018, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
 *
 * This file is part of HISAT 2.
 *
 * HISAT 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HISAT 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;

public class HT2Module {
	static {
		System.loadLibrary("ht2jni");
    }
    
    public static class HT2Position {
        int chr_id;
        int direction;
        long position;

        HT2Position(int id, int dir, int pos) {
            this.chr_id = id;
            this.direction = dir;
            this.position = pos;
        }

        @Override
        public String toString()
        {
            return getClass().getSimpleName() 
                    + "[chr_id=" + chr_id 
                    + ", direction=" + direction 
                    + ", position=" + position + "]";
        }
    }

    private native long init(String indexName, Map<String, Integer> options);
    private native void close(long handle);
    private native Map<String, Integer> get_options();

    private native String index_getrefnamebyid(long handle, int id);
    private native List<String> index_getrefnames(long handle);
    private native List<HT2Position> repeat_expand(long handle, String name, long rpos, long rlen);

    private long handle;

    HT2Module() {
        handle = 0;
    }

    public void initLibrary(String indexName) {
        initLibrary(indexName, null);
    }

    public void initLibrary(String indexName, Map<String, Integer> options)
    {
        if(handle == 0) {
            handle = init(indexName, options);
        }
    }

    public Map<String, Integer> initOption()
    {
        return get_options();
    }

    public String getRefNameById(int chr_id)
    {
        if(handle == 0) {
            // TODO: Exception
            return "";
        }

        // TODO: Exception(OutOfIndex)
        return index_getrefnamebyid(handle, chr_id);
    }

    public List<String> getRefNames()
    {
        if(handle == 0) {
            // TODO: Exception
            return new ArrayList<>();
        }
        return index_getrefnames(handle);
    }

    public List<HT2Position> repeatExpand(String name, long position, long length)
    {
        if(handle == 0) {
            // TODO: Exception
            return new ArrayList<>();
        }
        return repeat_expand(handle, name, position, length);
    }

    public void cleanup()
    {
        if(handle != 0) {
            close(handle);
            handle = 0;
        }
    }

    public void finalize() {
        cleanup();
    }

}
