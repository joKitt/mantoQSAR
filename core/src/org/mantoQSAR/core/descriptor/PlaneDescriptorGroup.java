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

import flexjson.JSON;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.mantoQSAR.core.DescriptorSet;
import org.mantoQSAR.core.Molecule;
import org.mantoQSAR.core.math.Plane;
import org.mantoQSAR.core.util.MoleculeTools;
import org.mantoQSAR.core.util.Projection;


public final class PlaneDescriptorGroup extends OrientationDescriptorGroup {


    public PlaneDescriptorGroup(){
        super(); 
        
        // number of descriptors in this group
        int nD = 32;

        this.descriptorList = new ArrayList<>();
        for (int i = 0; i < nD; i++) {
            this.descriptorList.add(new Descriptor("placeholder", 0.0));
        };

        
    }
    
    
    public PlaneDescriptorGroup(int observationID, int descriptorID) {
        super(observationID, descriptorID);
    }
    
    @JSON(include=false)
    private List<Double>getValueMap(){
        
        List<Double> valueMap = new ArrayList<>();
        DescriptorSet descriptorSet = this.getDescriptorSet();            
        Molecule m = this.getMolecule();
        List<Double[]> surfacePoints = m.getSurface();
        
        Projection p = new Projection();
        
        switch (descriptorSet.getSurface().getProperty()){
            case "esp":
                

                List<Double>  occ = m.getOccupancy(); 

                
                // assumes atom charges to be coded as occupancy value in PQR file 
                // according to PQR definition
                valueMap = p.calcProjection(surfacePoints, 
                        m.getAtomPosition(), 
                        m.getOccupancy(), 
                        descriptorSet.getSurface().getMapFunction(), 
                        0.0,            // minimum distance
                        78.0,           // dielectric constant of medium
                        this.getObservationSet().getCondition().getIonicStrength());
                break;
                
            case "hyd":
                valueMap = p.calcProjection(surfacePoints, 
                        m.getResiduePosition("mass"), 
                        m.getHydrophobicityConstant(), 
                        descriptorSet.getSurface().getMapFunction(), 0.0,            // minimum distance
                            78.0,           // dielectric constant of medium
                            this.getObservationSet().getCondition().getIonicStrength());
                break; 
                
            default: 
                logger.error("Property to be mapped in descriptor calculation not identified.");
                
                return valueMap;
        }
        
        p = null; 
        
        return valueMap; 
    }
    
    
    @Override
    public void calcDescriptor(){

        Molecule m = this.getMolecule();
        DescriptorSet descriptorSet = this.getDescriptorSet();            
        List<Double[]> surfacePoints = m.getSurface();
        
        logger.info("surface points identified " + surfacePoints.size());

                for(Double[] vector1:this.vector){
                    
                    List<Double[]> planeP = getPlane(vector1);
                    List<Double> valueMap2 = getPlaneValue(planeP, this.getValueMap());
                    Projection p = new Projection();
                    
                    Double[] surface2Plane = p.getAbsDistance(surfacePoints, planeP);
                    this.descriptorDetail.add((ArrayList<Descriptor>) this.calcDescriptorValue(valueMap2, surface2Plane));
                }
               
                logger.info("Select orientation based on " + descriptorSet.getProjection().getSelectID() + 
                 " selected for " + descriptorSet.getProjection().getSelectFunction());
                
                int descrInd = this.getDescriptorPosition(descriptorSet.getProjection().getSelectID());
                logger.info("Descriptor position in set is " + descrInd);
                
                int prefInd = this.getOrientation(descrInd, descriptorSet.getProjection().getSelectFunction());
                logger.info("Prefered orientation identified with " + prefInd);
                
                this.descriptorList = this.descriptorDetail.get(prefInd); 
    }

    private List<Descriptor> calcDescriptorValue(List<Double> valueMap, Double[] plane2surface) {

    List<Descriptor> descrList = new ArrayList<>();
    DescriptorSet descriptorSet = this.getDescriptorSet();
        
    String string = descriptorSet.getDescriptor().name;
    string = string.substring(0,1).toUpperCase() + string.substring(1);
    
    if(!string.substring(0,1).equalsIgnoreCase("_")){
        string = "_" + string;
    }

    List<Double> a_sasa = new ArrayList<>(); 
    int ii =0; 
        for(Double d:plane2surface){

            if (d < (10 + descriptorSet.getProjection().getDistance())) {
                a_sasa.add(valueMap.get(ii));
            }
            ii++;
        }
        double res = Math.pow(descriptorSet.getSurface().getResolution(), 2); 
    
        double sizeA = ((double) a_sasa.size())*(res*res);
        double size = (double) valueMap.size()*(res*res);

        descrList.add(new Descriptor("relSurfA" + string, sizeA/size));
        descrList.add(new Descriptor("totalSurfA" + string, sizeA));
        descrList.add(new Descriptor("nSurfP" + string, (double) a_sasa.size()));
        
        double sum = 0.0;
        for(Double dd:valueMap){
            sum = sum +dd; 
        }
        
        double sumA = 0.0; 
        for(Double dd:a_sasa){
            sumA = sumA+dd; 
        }
        
        double mean = sum/size;
        double meanA = sumA/sizeA;
        
        descrList.add(new Descriptor(("sum" + string), sum));
        descrList.add(new Descriptor(("sumSurfA" + string), sumA));
        descrList.add(new Descriptor(("densSurfA" + string), sumA/sizeA));
        descrList.add(new Descriptor(("mean" + string), mean));
        descrList.add(new Descriptor(("meanSurfA" + string), meanA));
        
        System.out.println( sizeA + "  " + meanA + ";...");
        
        
        List<Double> sortValueMap = new ArrayList<>(); 
        valueMap.stream().forEach((dd) -> {
            sortValueMap.add(dd);
        });
       
         List<Double> sortSurfA = new ArrayList<>(); 
         a_sasa.stream().forEach((dd) -> {
             sortSurfA.add(dd);
        });
        
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

    /*
    returns 3d coordinates of points in plane from selected index
    */
    @JSON(include=false)
    public List<Double[]> getPlane(int i){

        return getPlane(this.vector.get(i));
    }
    
    @JSON(include=false)
    public List<Double[]> getPlane(Double[] vector1){

        Molecule m = this.getMolecule(); 
        DescriptorSet descriptorSet = this.getDescriptorSet(); 
        
        Plane plane = new Plane(m.getAtomList(), vector1,
                           descriptorSet.getProjection().getSize(), 
                           descriptorSet.getProjection().getDensity(), 
                           descriptorSet.getProjection().getDistance());
        return plane.getPlane();
    }
    
    @JSON(include=false)
    public List<Double> getPlaneValue(int i){
        
        return getPlaneValue(getPlane(i), this.getValueMap()); 
    }
    
    @JSON(include=false)
    public List<Double> getPlaneValue(List<Double[]> plane, List<Double> valueMap){
       
         Molecule m = this.getMolecule();
         List<Double[]> surfacePoints = m.getSurface();
        
        Projection p = new Projection();
        
        return p.calcProjection(plane, 
                            surfacePoints, 
                            valueMap,
                            this.getDescriptorSet().getProjection().getMapFunction(),
                            0.0,            // minimum distance
                            78.0,           // dielectric constant of medium
                            this.getObservationSet().getCondition().getIonicStrength());
    }
    
    @JSON(include=false)
    public Double[] getCenterPoint() {
        return MoleculeTools.getCenter(this.getMolecule().getAtomList());
    } 
    
    @Override
    @JSON(include=false)
    public Integer getPreferredOrientationIndex() {
        DescriptorSet descriptorSet = this.getDescriptorSet();
        
        int descrInd = this.getDescriptorPosition(descriptorSet.getProjection().getSelectID());
                logger.info("Descriptor position in set is " + descrInd);
                
        return this.getOrientation(descrInd, descriptorSet.getProjection().getSelectFunction());
    }
}
