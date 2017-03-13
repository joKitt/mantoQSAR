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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.mantoQSAR.core.math.Vector;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class Projection {

    private static final Logger logger = LoggerFactory.getLogger(Projection.class); 
    Double ionicStrength; 
    Double dielectricConstant; 
    
    public Projection(){

    }
    
   
    public List<Double> calcProjection(List<Double[]> targetP, List<Double[]> sourceP, List<Double> value, int type, Double minDistance){
        
        Projection p = new Projection(); 
        return p.calcProjection(targetP, sourceP, value, type, minDistance, 2.0, 50.0); 
    }
    
    public List<Double> calcProjection(List<Double[]> targetP, 
                                              List<Double[]> sourceP, 
                                              List<Double> value, 
                                              int type, 
                                              Double minDistance, 
                                              Double dC, 
                                              Double iS) {

        this.ionicStrength = iS; 
        this.dielectricConstant = dC; 
        
        Double[] dist = new Double[3];
        List<Double> distM = new ArrayList<>();
        List<Double> mapValue = new ArrayList<>();
        
        if(targetP.isEmpty()){
            logger.error("Empty target list in Projection.calcProjection.");
            return null; 
        }
        
        if(sourceP.isEmpty()){
            logger.error("Empty source list in Projection.calcProjection.");
            return null; 
        }

        if (value.size() < sourceP.size()) {
            logger.error("Each source point should have a value. \n " + 
                    "sourceP " + sourceP.size() + 
                    "; valueP " + value.size());
            return mapValue;
        }
        
        for (int j = 0; j < targetP.size(); j++) {
          
            distM.clear();
            for (int k = 0; k < sourceP.size(); k++) {

                Double b = 0.0;
                for (int i = 0; i < dist.length; i++) {
                    Double a = (sourceP.get(k)[i] - targetP.get(j)[i]);
                    dist[i] = Math.pow(a, 2);
                    b = b + Math.pow(a, 2);
                }
                Double c = Math.sqrt(b);

                if (c < minDistance) {
                    c = minDistance;
                }

                distM.add(c);

            }

            // functions here 
            Double db = 0.304 / Math.sqrt(ionicStrength/1000); // debye screening length for monovalent ions
            Double v = 0.0;
            
            switch (type) {
                case 0:

                    for (int kk = 0; kk < distM.size(); kk++) {

                        v = v + (value.get(kk) * (1 / distM.get(kk)));
                    }
                    break;

                case 1:
                    for (int kk = 0; kk < distM.size(); kk++) {
                        v = v + (value.get(kk) * Math.pow(10, (-distM.get(kk) / 2)));
                    }
                    break;

                case 2:
                    for (int kk = 0; kk < distM.size(); kk++) {
                        v = v + (value.get(kk) * (1 / (1 + distM.get(kk))));
                    }
                    break;

                case 3:
                    for (int kk = 0; kk < distM.size(); kk++) {
                        v = v + (value.get(kk) * Math.pow(10, (-distM.get(kk))));
                    }
                    break;

                case 4:
                    for (int kk = 0; kk < distM.size(); kk++) {
                        v = v + (value.get(kk) * (1 / Math.pow(distM.get(kk), 6)));
                    }
                    break;

                case 5:

                    Double w = 1 / (4 * Math.PI * dielectricConstant);
                    for (int kk = 0; kk < distM.size(); kk++) {
                        v = v + (w * value.get(kk) * (1 / Math.pow(distM.get(kk), 2)));
                    }
                    break;

                case 6:

                    // distance given in Angstrom needs to be nm
                    for (int kk = 0; kk < distM.size(); kk++) {
                        v = v + ((value.get(kk)/ (dielectricConstant * (distM.get(kk)*0.1))) * Math.pow(10,((distM.get(kk)*-0.1) / db)) *1000) ;
                     }
                    break;
                    
                case 7:

                    Double wb = 1 / (4 * Math.PI * dielectricConstant);

                    for (int kk = 0; kk < distM.size(); kk++) {
                        v = v + (wb * value.get(kk) * (1 / Math.pow(distM.get(kk), 2)) * Math.pow(10, ((distM.get(kk)*-1) / db)) * 1000);
                    }
                    break;                  
            }
            mapValue.add(v);
        }
        return mapValue;
    }

    
    public Double[] getVectorDistance(List<Double[]> struct, Double[] vector, Double[] centerP) {
        // gives a distance value for each point in struct list, distance measured along a vector; 

        Double[] d2Point = new Double[struct.size()];

        for (int i = 0; i < struct.size(); i++) {

            Double[] nSurfPVec = new Double[3];
            nSurfPVec[0] = struct.get(i)[0] - centerP[0];
            nSurfPVec[1] = struct.get(i)[1] - centerP[1];
            nSurfPVec[2] = struct.get(i)[2] - centerP[2];

            Double c = 0.0;
            nSurfPVec[0] = nSurfPVec[0] / Vector.norm(vector);
            c = c + Math.abs(nSurfPVec[0] - vector[0]);
            nSurfPVec[1] = nSurfPVec[1] / Vector.norm(vector);
            c = c + Math.abs(nSurfPVec[1] - vector[1]);
            nSurfPVec[2] = nSurfPVec[2] / Vector.norm(vector);
            c = c + Math.abs(nSurfPVec[2] - vector[2]);

            d2Point[i] = c;
        }
        return d2Point;
    }
    
    /*
    * absolut distance of object points to plane
    */
    public Double[] getAbsDistance(List<Double[]> planeP, List<Double[]> struct) {
        
        Double[] d2Plane = new Double[struct.size()];
        List<Double> DistA = new ArrayList<>();
        for (int i = 0; i < struct.size(); i++) {

            DistA.clear();
            
            for (Double[] planePoint : planeP) {
                Double b = Math.sqrt(Math.pow(struct.get(i)[0] - planePoint[0], 2) 
                                   + Math.pow(struct.get(i)[1] - planePoint[1], 2) 
                                   + Math.pow(struct.get(i)[2] - planePoint[2], 2));
                DistA.add(b);
            }
            
            if(DistA.size() > 1){
            d2Plane[i] = Collections.min(DistA);
            }
            if(DistA.size()==1){
                d2Plane[i] = DistA.get(0);
            }
        }
        return d2Plane;
    }
    
     public static Double getMinAbsDistance(List<Double[]> planeP, List<Double[]> struct) {

        Projection project = new Projection(); 
        Double[] d2Plane = project.getAbsDistance(planeP, struct);
        Double dist = Collections.min(Arrays.asList(d2Plane));
        return dist;
    }
}
