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
import org.mantoQSAR.core.util.Projection;

public final class SurfaceDescriptorGroup extends DescriptorGroup{

    List<Double> valueMap; 

    public SurfaceDescriptorGroup(){
        super();
    }
    
    public SurfaceDescriptorGroup(int observationID, int descriptorID) {
        super(observationID, descriptorID);
        
         this.valueMap = new ArrayList<>();
   
    // number of descriptors in this group
        int nD = 28;

        this.descriptorList = new ArrayList<>();
        for (int i = 0; i < nD; i++) {
            this.descriptorList.add(new Descriptor("placeholder", 0.0));
        };
    }

    @Override
    public void calcDescriptor(){
        
        Molecule m = this.getMolecule();
        DescriptorSet descriptorSet = this.getDescriptorSet();             
        List<Double[]> surfacePoints = m.getSurface();

        this.valueMap = new ArrayList<>(surfacePoints.size());
        
        logger.info("surface points identified " + surfacePoints.size());
        Projection p = new Projection(); 
        
        switch (descriptorSet.getSurface().getProperty()){
            case "esp":
                
                // assumes atom charges to be coded as occupancy value in PQR file 
                // according to PQR definition
                this.valueMap = p.calcProjection(surfacePoints, 
                        m.getAtomPosition(), 
                        m.getOccupancy(), 
                        descriptorSet.getSurface().getMapFunction(), 0.0);
                break;
                
            case "hyd":
                this.valueMap = p.calcProjection(surfacePoints, 
                        m.getResiduePosition("mass"), 
                        m.getHydrophobicityConstant(), 
                        descriptorSet.getSurface().getMapFunction(), 0.0);
                break; 
                
            default: 
                logger.error("Property to be mapped in descriptor calculation not identified.");
                return;
        }
        
        p = null; 
                this.descriptorList = this.calcDescriptorValue(valueMap);
                
    }

    
private List<Descriptor> calcDescriptorValue(List<Double> valueMap) {
      
    List<Descriptor> descrList = new ArrayList<>();
    DescriptorSet descriptorSet = this.getDescriptorSet();
        
    stIdent = this.stIdent.substring(0,1).toUpperCase() + this.stIdent.substring(1); 
    stIdent = "_" + stIdent; 
    
        double size = (double) valueMap.size();
        double res = Math.pow(descriptorSet.getSurface().getResolution(), 2); 

        double c = size / res; 
        descrList.add(new Descriptor("totalSurf", c));
        descrList.add(new Descriptor("nSurfP", size));
        
        double sum = 0.0;
        for(Double dd:valueMap){
            sum = sum +dd; 
        }
        
        double mean = sum/size;
        descrList.add(new Descriptor(("sum" + stIdent), sum));
        descrList.add(new Descriptor(("mean" + stIdent), mean));
    
        List<Double> sortValueMap = new ArrayList<>(); 
        for (Double dd:valueMap){
            sortValueMap.add(dd); 
        }

        Collections.sort(sortValueMap); 
       
        int center = sortValueMap.size()/2; 
        double median = sortValueMap.get(center);
        double max = Collections.max(valueMap);
        double min = Collections.min(valueMap);
        
        descrList.add(new Descriptor(("median" + stIdent), median));
        descrList.add(new Descriptor(("meanRes" + stIdent), (mean)/res));
        
        descrList.add(new Descriptor(("max" + stIdent),max));
        descrList.add(new Descriptor(("min" + stIdent),min));
        
        descrList.add(new Descriptor(("devA" + stIdent),(max - min)/median));
        descrList.add(new Descriptor(("devB" + stIdent),(max - min)/mean));
        
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

        descrList.add(new Descriptor(("nPos" + stIdent), (double) posValueMap.size()));
        descrList.add(new Descriptor(("nNeg" + stIdent), (double) posValueMap.size()));
        
        descrList.add(new Descriptor(("relPos" + stIdent), (double) posValueMap.size()/ valueMap.size()));
        descrList.add(new Descriptor(("relNeg" + stIdent), (double) posValueMap.size()/ valueMap.size()));
        
        descrList.add(new Descriptor(("sumPos" + stIdent), sumPosVal));
        descrList.add(new Descriptor(("sumNeg" + stIdent), sumNegVal));
        
        descrList.add(new Descriptor(("averPos" + stIdent), sumPosVal/ posValueMap.size()));
        descrList.add(new Descriptor(("averNeg" + stIdent), sumNegVal/ negValueMap.size()));

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
            descrList.add(new Descriptor(("binAbs" + stIdent + "_" + i), (double) binVal.get(i)/res));
        }

        return descrList; 
} 

public List<Double> getSurfaceValue(){
    
    return this.valueMap;
}

public void setSurfaceValue(List<Double> vm){
    
    this.valueMap = vm; 
}
}
