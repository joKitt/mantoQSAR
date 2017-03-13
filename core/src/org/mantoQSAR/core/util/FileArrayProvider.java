/* This file is part of mantoQSAR.

 mantoQSAR - Quantitative structure-activity relationship descriptor 
 calculation and modeling for biomolecules.
			
 Copyright (C) 2016  Jörg Kittelmann


 mantoQSAR is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, 
 or any later version.

 mantoQSAR is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with mantoQSAR. If not, see <http://www.gnu.org/licenses/>.
 */
package org.mantoQSAR.core.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class FileArrayProvider {

    public List<Double[]> readLines(String filename) throws IOException {
        FileReader fileReader = new FileReader(filename);

        List<String> lines;
        String line = null;

        try (BufferedReader bufferedReader = new BufferedReader(fileReader)) {
            lines = new ArrayList<String>();

            while ((line = bufferedReader.readLine()) != null) {
                // remove parenthesis 
                line = line.replace("{", "");
                line = line.replace("}", "");

                lines.add(line);
            }
        }
        return toDouble(lines);
    }

    public List<Double[]> toDouble(List<String> lines) {

        List<Double[]> valList = new ArrayList<Double[]>();

        for (String string : lines) {
            // split string 
            String[] s = string.split(" ");
            Double[] val = new Double[3];

            for (int i = 0; i < s.length; i++) {
                val[i] = Double.parseDouble(s[i]);
            }
            valList.add(val);
        }

        return valList;
    }
}
