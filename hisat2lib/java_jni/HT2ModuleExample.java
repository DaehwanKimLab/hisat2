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

public class HT2ModuleExample {

	public static void main(String[] args) {
        HT2Module module = null;
        Map<String, Integer> ht2Options = null;
        String indexPath = "../../evaluation/indexes/HISAT2_22/22_rep";
        //String indexPath = "../../evaluation/indexes/HISAT2/genome_rep";

        try {
            module = new HT2Module();            
            
            // Get Default Options
            //ht2Options = module.InitOption();
            // or
            ht2Options = new HashMap<>();

            //ht2Options.put("gVerbose", 1);
            //ht2Options.put("startVerbose", 1);
            //ht2Options.put("sanityCheck", 1);

            System.out.println(ht2Options);

            module.InitLibrary(indexPath, ht2Options);
            //module.InitLibrary(indexPath);

            System.out.println("CHR_ID 0: " + module.GetRefNameById(0));

            // in 22_rep, OutOfIndex
            System.out.println("CHR_ID 1: " + module.GetRefNameById(1));

            List<String> refnames = module.GetRefNames();

            System.out.println("Refnames size:" + refnames.size());
            for(String name: refnames) {
                System.out.println(name);
            }
           
            // repeat for 22_rep
            List<HT2Module.HT2Position> positions = module.RepeatExpand("rep100-300", 8308, 100);
            
            // repeat for genome_rep
            //List<HT2Module.HT2Position> positions = module.RepeatExpand("rep100-300", 2446692, 100);

            System.out.println("Repeat expand size: " + positions.size());
            for(HT2Module.HT2Position pos: positions) {
                String chrName = refnames.get(pos.chr_id).split(" ")[0];
                String direction = pos.direction == 0 ? "+":"-";

                System.out.println(chrName + ":" + pos.position + ":" + direction); 
            }

        } finally {
            if(module != null) {
                module.Cleanup();
            }
        }
	}
}
