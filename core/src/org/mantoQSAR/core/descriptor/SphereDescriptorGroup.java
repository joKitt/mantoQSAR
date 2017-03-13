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


package org.mantoQSAR.core.descriptor;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.mantoQSAR.core.DescriptorSet;
import org.mantoQSAR.core.Molecule;
import org.mantoQSAR.core.math.Plane;
import org.mantoQSAR.core.util.Projection;

public final class SphereDescriptorGroup extends OrientationDescriptorGroup {

    Double ionicStrength = 0.0;
    Double dielectricConstant = 2.0;
    Double minDistance = 0.0;
   
    public SphereDescriptorGroup(){
        super();
    }
    
    public SphereDescriptorGroup(int observationID, int descriptorID) {
        super(observationID, descriptorID);

        // number of descriptors in this group
        int nD = 34;

        this.descriptorList = new ArrayList<>();
        for (int i = 0; i < nD; i++) {
            this.descriptorList.add(new Descriptor("placeholder", 0.0));
        };      
    }

    @Override
    public void calcDescriptor(){

        Molecule m = this.getMolecule();
        DescriptorSet descriptorSet = this.getDescriptorSet();
        
        List<Double> valueMap;              
        List<Double[]> surfacePoints = m.getSurface();

        logger.info("surface points identified " + surfacePoints.size());
        
        Projection p = new Projection(); 
        
        switch (descriptorSet.getSurface().getProperty()){
            case "esp":
                
                // assumes atom charges to be coded as occupancy value in PQR file 
                // according to PQR definition
                valueMap = p.calcProjection(surfacePoints, 
                        m.getAtomPosition(), 
                        m.getOccupancy(), 
                        descriptorSet.getSurface().getMapFunction(), 0.0);
                 
                break;
                
            case "hyd":
                valueMap = p.calcProjection(surfacePoints, 
                        m.getResiduePosition("mass"), 
                        m.getHydrophobicityConstant(), 
                        descriptorSet.getSurface().getMapFunction(), 0.0);
                break; 
                
            default: 
                logger.error("Property to be mapped in descriptor calculation not identified.");
                
                return;
        }
        
        p = null;

                for(Double[] vector1:this.vector){
                   Plane plane = new Plane(m.getAtomList(), vector1,
                           descriptorSet.getProjection().getSize(), 
                           descriptorSet.getProjection().getDensity(), 
                           descriptorSet.getProjection().getDistance());
                    List<Double[]> planeP = plane.getPlane();
                
                    logger.info("plane points identified " + planeP.size());
                    
                    p = new Projection(); 
                    
                    List<Double> valueMap2 = p.calcProjection(planeP, 
                            surfacePoints, 
                            valueMap, 
                            descriptorSet.getProjection().getMapFunction(), 
                            0.0,            // minimum distance
                            78.0,           // dielectric constant of medium
                            this.getObservationSet().getCondition().getIonicStrength());
                    
                    p = null; 
                    logger.info("plane values mapped " + valueMap2.size());
                    
                    p = new Projection();
                    Double[] surface2Plane = p.getAbsDistance(surfacePoints, planeP);
                    this.descriptorDetail.add((ArrayList<Descriptor>) this.calcDescriptorValue(valueMap2, surface2Plane));
                }
               
                this.descriptorList = this.calcAverage(); 
               
    }

    private List<Descriptor> calcDescriptorValue(List<Double> valueMap, Double[] plane2surface) {

        DescriptorSet descriptorSet = this.getDescriptorSet();
        
    List<Descriptor> descrList = new ArrayList<>();
        
    String string = descriptorSet.getDescriptor().name;
    string = string.substring(0,1).toUpperCase() + string.substring(1);
    
    if(!string.substring(0,1).equalsIgnoreCase("_")){
        string = "_" + string;
    }

    List<Double> a_sasa = new ArrayList<>(); 
        for(Double d:plane2surface){
            if (d < (2.8 + descriptorSet.getProjection().getDistance())) {
                a_sasa.add(d);
            }
        }
    
        double sizeA = (double) a_sasa.size();
        double size = (double) valueMap.size();
        double res = Math.pow(descriptorSet.getSurface().getResolution(), 2); 
        
        
        descrList.add(new Descriptor("nSurfA" + string, sizeA));
 
        descrList.add(new Descriptor("relSurfA" + string, sizeA/size));
        
        double c = sizeA / res; 
        descrList.add(new Descriptor("totalSurfA" + string, c));
        
        descrList.add(new Descriptor("nSurfP" + string, size));
        
        double sum = 0.0;
        for(Double dd:valueMap){
            sum = sum +dd; 
        }
        
        double sumA = 0.0; 
        for(Double dd:a_sasa){
            sumA = sumA +dd; 
        }
        
        
        double mean = sum/size;
        double meanA = sumA/sizeA;
        
        descrList.add(new Descriptor(("sum" + string), sum));
        descrList.add(new Descriptor(("sumSurfA" + string), sumA));
        
        descrList.add(new Descriptor(("mean" + string), mean));
        descrList.add(new Descriptor(("meanA" + string), meanA));
        
        List<Double> sortValueMap = new ArrayList<>(); 
        for (Double dd:valueMap){
            sortValueMap.add(dd); 
        }
       
        List<Double> sortSurfA = new ArrayList<>(); 
        for (Double dd:a_sasa){
            sortSurfA.add(dd); 
        }
        
        Collections.sort(sortValueMap); 
        Collections.sort(sortSurfA); 
        
        int center = sortValueMap.size()/2; 
        int centerA = sortSurfA.size()/2;
        
        double median = sortValueMap.get(center);
        double medianA = sortValueMap.get(centerA);
        
        double max = Collections.max(valueMap);
        double min = Collections.min(valueMap);
        
        descrList.add(new Descriptor(("median" + string), median));
        descrList.add(new Descriptor(("medianSurfA" + string), medianA));
        descrList.add(new Descriptor(("meanRes" + string), (mean)/res));
        descrList.add(new Descriptor(("meanResA" + string), (meanA)/res));
        
        descrList.add(new Descriptor(("max" + string),max));
        descrList.add(new Descriptor(("min" + string),min));
        
        descrList.add(new Descriptor(("devA" + string),(max - min)/median));
        descrList.add(new Descriptor(("devB" + string),Collections.max(valueMap)));
        
        List<Double> posValueMap = new ArrayList<>(); 
        List<Double> negValueMap = new ArrayList<>();
        double sumPosVal = 0.0;
        double sumNegVal = 0.0;
        
        for (Double dd:valueMap){
            if(dd > 0.0){
            posValueMap.add(dd);
            sumPosVal = sumPosVal + dd;}
            if(dd < 0.0){
            negValueMap.add(dd);
            sumNegVal=sumNegVal + dd;
            }
        }
        
        descrList.add(new Descriptor(("nPos" + string), (double) posValueMap.size()));
        descrList.add(new Descriptor(("nNeg" + string), (double) posValueMap.size()));
        
        descrList.add(new Descriptor(("relPos" + string), (double) posValueMap.size()/ valueMap.size()));
        descrList.add(new Descriptor(("relNeg" + string), (double) posValueMap.size()/ valueMap.size()));
        
        descrList.add(new Descriptor(("sumPos" + string), sumPosVal));
        descrList.add(new Descriptor(("sumNeg" + string), sumNegVal));
        
        descrList.add(new Descriptor(("averPos" + string), sumPosVal/ posValueMap.size()));
        descrList.add(new Descriptor(("averNeg" + string), sumNegVal/ negValueMap.size()));
        
        double binScale = descriptorSet.getDescriptor().binScale;
        int nBin = 10; 
        List<Double> binVal = new ArrayList<>();
   
   
        for(int i =1; i <= nBin; i++){
            double binLow= (-(binScale*nBin)*0.5) + (i)* binScale;
            double binHigh = (-(binScale*nBin)*0.5) + (i+1)* binScale;
            
            double a = 0.0;
        for (Double valueMap1 : valueMap) {
            if (valueMap1 >= binLow && valueMap1 < binHigh) {
                a = a+1;   
            }
        }
        binVal.add(a);
        }
        
        for(int i = 0; i < nBin; i++){ 
            descrList.add(new Descriptor(("binAbs" + string + "_" + i), (double) binVal.get(i)));
        }
        return descrList; 
    }

    private List<Descriptor> calcAverage() {
        
        List<Descriptor> descrList = new ArrayList<>();
        int m = this.descriptorDetail.size();
        int n = this.descriptorDetail.get(0).size();
        
        for(int i = 0; i < n; i++){
            Descriptor d = new Descriptor(this.descriptorDetail.get(0).get(i).getName(), 0.0);
            double dd = 0.0; 

            for(int j = 0; j < m; j++){
                try{
                dd += this.descriptorDetail.get(j).get(i).getValue();
                }catch(Exception e){
                    logger.error(e.getMessage());
                }
            }

            d.value = dd/m; 
            descrList.add(d);
        }
        
        return descrList;
    }

    public Double getIonicStrength() {
        return ionicStrength;
    }

    public void setIonicStrength(Double ionicStrength) {
        this.ionicStrength = ionicStrength;
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
