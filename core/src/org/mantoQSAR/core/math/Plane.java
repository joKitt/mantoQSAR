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

package org.mantoQSAR.core.math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;
import org.biojava.nbio.structure.Atom;
import org.mantoQSAR.core.util.MoleculeTools;
import org.mantoQSAR.core.util.Projection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Plane{
    
   Logger logger;
   List<Double[]> planePoints;  
   Double dielectricConstant = 2.0; 
   Double minDistance = 0.0; 
   String type; 
   
   Double[] vector; 
   Double size = 100.0; 
   Double density = 5.0; 
   Double distance;
   
   boolean atom = false;
   
   
   List<Atom> atomList;

    public Plane(List<Atom> atomList, Double[] vector, Double size, Double density, Double distance){
        
        logger = LoggerFactory.getLogger(Plane.class);
        
        atom = true; 
        this.atomList = atomList;
        this.vector = vector;
        this.size = size; 
        this.density = density; 
        this.distance = distance; 
        planePoints = new ArrayList<>(); 

    }

    public List<Double[]> getPlane(){
        planePoints.clear();
        planePoints = this.calcPlane(this.atomList, this.vector, this.size, density, distance);
        return planePoints; 
    }
    
    public List<Double[]> calcPlane(List<Atom> atomList, Double[] vector, Double size, Double density, Double distance) {

        List<Double[]> planeP = new ArrayList<>();
         
        Double oldSize = size;
        Double dp = Math.ceil(size / density);
        size = density * dp;

        if (!oldSize.equals(size)) {
            logger.info("Plane size recalculated to " + size.toString());
        }
        
        
        List<Double[]> refP = MoleculeTools.getAtomPosition(atomList);

        Double[] cP = MoleculeTools.getCenter(atomList);

        for (int i = 0; i < vector.length; i++) {
            vector[i] = vector[i] / Vector.norm(vector);
        }
            
       // center of plane to be build
        Double[] pSn = new Double[3];

        for (int i = 0; i < 3; i++) {
            pSn[i] = cP[i] + (vector[i] * 400);
        }
        // use start vectors in 3 directions to reduce error by limitation of double number
        Double[] vec1a = new Double[]{1.0, 1.0, 0.0};
        Double[] vec1b = new Double[]{1.0, 0.0, 0.0};
        Double[] vec1c = new Double[]{0.0, 1.0, 0.0};
        vec1a[2] = (((vec1a[0] * vector[0]) + (vec1a[1] * vector[1])) / vector[2])*-1;
        vec1b[2] = (((vec1b[0] * vector[0]) + (vec1b[1] * vector[1])) / vector[2])*-1;
        vec1c[2] = (((vec1c[0] * vector[0]) + (vec1c[1] * vector[1])) / vector[2])*-1;
        
        // norm vec1s before use
        for (int i = 0; i < vec1a.length; i++) {
            vec1a[i] = vec1a[i] / Vector.norm(vec1a);
            vec1b[i] = vec1b[i] / Vector.norm(vec1b);
            vec1c[i] = vec1c[i] / Vector.norm(vec1c);
        }

        Double[] vec2a = Vector.crossProduct(vector, vec1a);
        Double[] vec2b = Vector.crossProduct(vector, vec1b);
        Double[] vec2c = Vector.crossProduct(vector, vec1c);
        
        Double[] vec3a = Vector.crossProduct(vec1a, vec2a);
        Double[] vec3b = Vector.crossProduct(vec1b, vec2b);
        Double[] vec3c = Vector.crossProduct(vec1c, vec2c);
        
         // norm vec1s before use
        for (int i = 0; i < vec3a.length; i++) {
            vec3a[i] = vec3a[i] / Vector.norm(vec3a);
            vec3b[i] = vec3b[i] / Vector.norm(vec3b);
            vec3c[i] = vec3c[i] / Vector.norm(vec3c);
        }
        
        // difference between start vector and vector of back projection
        Double aa = Math.abs(vec3a[0]-vector[0]) + Math.abs(vec3a[1]-vector[1]) + Math.abs(vec3a[2]-vector[2]);
        Double ab = Math.abs(vec3b[0]-vector[0]) + Math.abs(vec3b[1]-vector[1]) + Math.abs(vec3b[2]-vector[2]);
        Double ac = Math.abs(vec3c[0]-vector[0]) + Math.abs(vec3c[1]-vector[1]) + Math.abs(vec3c[2]-vector[2]);

        Double[] vec1;
        Double[] vec2;
        
        
        if (aa <= ab & aa <= ac) {
            // take vector a
            vec1 = ArrayUtils.clone(vec1a);
            vec2 = ArrayUtils.clone(vec2a);
        } else {
            if (ab <= aa & ab <= ac) {
                // take vector b
                vec1 = ArrayUtils.clone(vec1b);
                vec2 = ArrayUtils.clone(vec2b);
            } else {
                // take vector c
                vec1 = ArrayUtils.clone(vec1c);
                vec2 = ArrayUtils.clone(vec2c);
            }
        }

        // for debuging cross product of vec1 and vec2 should be vector
        Double[] vecCross= Vector.crossProduct(vec1a, vec2);
        // norm vec1 before use
        for (int i = 0; i < vecCross.length; i++) {
            vecCross[i] = vecCross[i] / Vector.norm(vecCross);
        }
        
        for (int i = 0; i < vec1.length; i++) {
            vec1[i] = vec1[i] / Vector.norm(vec1);
            vec2[i] = vec2[i] / Vector.norm(vec2);

        }
        
        Double[] pPlaneEdge = new Double[3];

        // get edge point of plane
        for (int i = 0; i < pSn.length; i++) {
            pPlaneEdge[i] = pSn[i] - (vec1[i] * (size/2)) - (vec2[i] * (size/2)); 
        }
        Double count = (size / density); 
        List<Double[]> M = new ArrayList<>(); 

        for(int i = 0; i < count; i++){
            for(int j = 0; j < count; j++){
                
                Double[] point = new Double[3];
                for(int k =0; k < 3; k++){
                    point[k] = pPlaneEdge[k] + vec1[k]*density*i + vec2[k]*density*j; 
                }
                M.add(point);
             }
        }
 
        Double minD2Plane = Projection.getMinAbsDistance(M, refP);

        for (int i = 0; i < vector.length; i++) {
            vector[i] = vector[i] / Vector.norm(vector);
        }
        
        Double[] corVec = new Double[]{vector[0] * minD2Plane, vector[1] * minD2Plane, vector[2] * minD2Plane};
       
        Double corLength = Math.sqrt(Math.pow(corVec[0], 2) + Math.pow(corVec[1], 2) + Math.pow(corVec[2], 2));
       
        List<Double[]> planeCorrected = new ArrayList<>();
        
        for (Double[] M1 : M) {
            Double[] p = new Double[3];
            p[0] = M1[0] - corVec[0];
            p[1] = M1[1] - corVec[1];
            p[2] = M1[2] - corVec[2];
            planeCorrected.add(p);
        }

        // add distance according to parameter set
        Double corDistance = distance - (minD2Plane - corLength); 
        Double[] distVec = new Double[]{vector[0] * corDistance, vector[1] * corDistance, vector[2] * corDistance};

        planeP.clear();
        for (Double[] M1 : planeCorrected) {
            Double[] pp = new Double[3];
            pp[0] = M1[0] + distVec[0];
            pp[1] = M1[1] + distVec[1];
            pp[2] = M1[2] + distVec[2];
            planeP.add(pp);
        }
        
        return planeP; 
    }

    public Double[] getAbsDistance(List<Double[]> struct){ 
        Projection p = new Projection();
        return p.getAbsDistance(this.planePoints, struct); 
    }
  
    public Double getMinAbsDistance(List<Double[]> struct) {

        Double[] d2Plane = this.getAbsDistance(struct);
        Double dist = Collections.min(Arrays.asList(d2Plane));
        return dist;
    }

    public Double getDielectricConstant() {
        return dielectricConstant;
    }

    public void setDielectricConstant(Double dielectricConstant) {
        this.dielectricConstant = dielectricConstant;
    }

    public Double getMinDistance() {
        return minDistance;
    }

    public void setMinDistance(Double minDistance) {
        this.minDistance = minDistance;
    }
    
}
