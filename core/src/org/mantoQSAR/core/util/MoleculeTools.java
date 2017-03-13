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

import java.util.ArrayList;
import java.util.List;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;



public abstract class MoleculeTools {

    private static final Logger logger = LoggerFactory.getLogger(MoleculeTools.class);
    
    public static Double[] getCenter(List<Atom> aList) {

        List<Double> weight = new ArrayList<>();
        for (Atom aList1 : aList) {
            weight.add(1.0);
        }
        return MoleculeTools.getCenter(aList, weight);
    }

    public static Double[] getCenter(List<Atom> aList, List<Double> weight) {

        Double[] center = new Double[3];

        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double w = 0.0;

        for (int i = 0; i < aList.size(); i++) {

            x = x + aList.get(i).getX() * weight.get(i);
            y = y + aList.get(i).getY() * weight.get(i);
            z = z + aList.get(i).getZ() * weight.get(i);
            w = w + weight.get(i);
        }

        Double wAv = w/aList.size(); 
        center[0] = x / (aList.size() * wAv);
        center[1] = y / (aList.size() * wAv);
        center[2] = z / (aList.size() * wAv);

        return center;
    }

    public static List<Double> getAtomMass(List<Atom> aList) {
        List<Double> atomMass = new ArrayList<>();

        for (Atom a : aList) {
            atomMass.add(new Double(a.getElement().getAtomicMass()));
        }
        return atomMass;
    }

    public static List<Double> getHydrophobicityConstant(List<Group> gList) {

        List<Double> hydConst = new ArrayList<>();

        for (Group g : gList) {
          
            String code = g.getChemComp().getThree_letter_code();
            
            logger.trace("hydrophobicity constant for " + g.getChemComp().getThree_letter_code());
            /* values according to:
             Kyte, J., & Doolittle, R. F. (1982). A simple method  for displaying the
             hydropathic character of a protein.  Journal of molecular biology, 157(1), 105-32.
             http://www.ncbi.nlm.nih.gov/pubmed/7108955 */
            
            Double val;
            
            if(code != null){
                
            switch (code) {
                case "TRP":
                    val = 7.9;
                    break;
                case "PHE":
                    val = 7.5;
                    break;
                case "LEU":
                    val = 6.6;
                    break;
                case "ILE":
                    val = 4.3;
                    break;
                case "TYR":
                    val = 7.1;
                    break;
                case "VAL":
                    val = 5.1;
                    break;
                case "MET":
                    val = 2.5;
                    break;
                case "PRO":
                    val = 2.2;
                    break;
                case "CYS":
                    val = 0.0;
                    break;
                case "ARG":
                    val = -1.1;
                    break;
                case "ALA":
                    val = -0.3;
                    break;
                case "LYS":
                    val = -3.6;
                    break;
                case "GLY":
                    val = 1.2;
                    break;
                case "ASP":
                    val = -1.4;
                    break;
                case "GLU":
                    val = 0.0;
                    break;
                case "HIS":
                    val = -1.3;
                    break;
                case "THR":
                    val = -2.2;
                    break;
                case "SER":
                    val = -0.6;
                    break;
                case "ASN":
                    val = -0.2;
                    break;
                case "GLN":
                    val = -0.2;
                    break;
                default:
                    val = 0.0;
                    break;
                }
            }else{
                val = 0.0; 
                }
            hydConst.add(val);
        }

        return hydConst;
    }
    
    public static List<Double[]> getAtomPosition(List<Atom> aList) {

        List<Double[]> atomPos = new ArrayList<>(aList.size());

        for (Atom a : aList) {
            double[] aP = a.getCoords();
            Double[] aP1 = new Double[3];
            for (int i = 0; i < aP.length; i++) {
                aP1[i] = (Double) aP[i];
            }
            atomPos.add(aP1);
        }
        return atomPos;
    }
}
