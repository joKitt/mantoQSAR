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


package org.mantoQSAR.core;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Reader;
import java.lang.reflect.Type;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.mantoQSAR.core.descriptor.*;
import org.mantoQSAR.core.math.Matrix;
import org.mantoQSAR.core.math.Vector;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.mantoQSAR.core.descriptor.Descriptor;
import org.mantoQSAR.core.util.ColorStatic;


public final class Screen {

    Logger logger;
    private static Screen instance = null;
    final private PropertyChangeSupport changes;

    Matrix modelDescriptorMatrix = new Matrix(3, 20);
    Matrix predictDescriptorMatrix = new Matrix(3, 20);
    Matrix descriptorMatrix = new Matrix(3,3);

    List<ObservationSet> observationSetList;
    List<DescriptorSet> descriptorSetList;
    List<Molecule> moleculeList;
    ProjectDescriptor projectDescriptor; 
    
    Boolean idle = true;
    
    // fixed value for development purpose
    Integer settingProcessor = 0; // 1 = multi thread; 0 = single thread
    String projectPath;

    public static Screen getInstance() {

        if (instance == null) {
            instance = new Screen();
        }
        return instance;
    }

    protected Screen() {
        this.changes = new PropertyChangeSupport(this);

        logger = LoggerFactory.getLogger(Screen.class);

        descriptorSetList = new ArrayList<>();
        observationSetList = new ArrayList<>();
        moleculeList = new ArrayList<>();
        this.projectDescriptor = new ProjectDescriptor();

        projectPath = System.getProperty("user.home") + "/mantoQSARProjects";
    }

    public void loadDescriptorSetting(final String s) {

        logger.info("Reading descriptor settings from file " + s);
        System.out.println(ColorStatic.PURPLE + "Reading descriptor settings from file " + s + ColorStatic.RESET);

        File sFile = new File(s);
        if (!sFile.exists()) {
            logger.error("Descriptor setting file " + sFile.getAbsolutePath() + " not found.");
            return;
        }

        Thread t1 = new Thread(new Runnable() {

            @Override
            public void run() {

                try {
                    Reader reader;
                    reader = new BufferedReader(new FileReader(new File(s)));

                    Gson gson = new GsonBuilder().create();
                    Type listType = new TypeToken<ArrayList<DescriptorSet>>() {
                    }.getType();

                    descriptorSetList = gson.fromJson(reader, listType);

                } catch (FileNotFoundException ex) {
                    logger.error(ex.getMessage());
                }
                changes.firePropertyChange(MoleculeStatic.CHANGE_LOAD_DESCRIPTOR_SETTING, null, null);

            }
        });
        t1.start();

        descriptorSetList.stream().forEach((dS) -> {
            logger.debug(dS.toString());
        });
    }

    private void displayDescriptorSetting() {
        for (DescriptorSet ds : this.descriptorSetList) {
            System.out.println("name \t\t" + ds.name + "\n"
                    + "surface \n"
                    + "   resolution \t" + Float.toString(ds.surface.resolution) + "\n"
                    + "   type of surface \t" + ds.getSurface().getTypeOfSurface() + "\n"
                    + "   size of sphere \t" + ds.getSurface().getSizeOfSphere() + "\n"
                    + "   property \t\t" + ds.surface.getProperty() + "\n"
                    + "   map function \t" + Integer.toString(ds.surface.mapfunc));
            if (ds.projection != null) {
                System.out.println("projection \n"
                        + "   type \t\t" + ds.projection.getType() + "\n"
                        + "   density \t\t" + ds.projection.density.toString() + "\n"
                        + "   size \t\t" + ds.projection.size.toString() + "\n"
                        + "   distance \t\t" + ds.projection.distance.toString() + "\n"
                        + "   map function \t" + Integer.toString(ds.projection.mapfunc) + "\n"
                        + "   orientations \t" + Integer.toString(ds.projection.orientation) + "\n");
            }
            System.out.println("descriptor \n"
                    + "   name \t\t" + ds.descriptor.name);

            if (ds.descriptor.property != null) {
                System.out.println("   property \t\t" + ds.descriptor.property);
            }
            if (ds.descriptor.binScale != null) {
                System.out.println("   binScale \t\t" + ds.descriptor.binScale.toString());
            }
        }

    }

    public void loadObservationSetting(final String s) {

        logger.info("Reading observation settings from file " + s);
        System.out.println(ColorStatic.PURPLE + "Reading observation settings from file " + s + ColorStatic.RESET);

        File sFile = new File(s);

        if (!sFile.exists()) {
            logger.error("Observation setting file " + sFile.getAbsolutePath() + " not found.");
            return;
        }

        Thread t1 = new Thread(new Runnable() {
            @Override
            public void run() {

                try {
                    observationSetList.clear();
                    moleculeList.clear();

                    ImportExport io = new ImportExport();  

                    Reader reader;
                    reader = new BufferedReader(new FileReader(new File(s)));

                    Gson gson = new GsonBuilder().create();
                    Type listType = new TypeToken<ArrayList<ObservationSet>>() {
                    }.getType();
                    observationSetList = gson.fromJson(reader, listType);

                    // import structures from PDB/PQR files
                    for (ObservationSet obs : observationSetList) {
                       logger.info("Import structure for " + obs.getName());

                        Molecule m = new Molecule();
                        System.out.println("reading structure from file " + obs.getFile());
                        
                        
                        m.setStructure(io.importStructure(new File(projectPath + obs.getFile()))); 

                        // center molecule in coord system
                        m.setIonicStrength(obs.getCondition().ionicStrength);
                        moleculeList.add(m);

                    }
                    changes.firePropertyChange(MoleculeStatic.CHANGE_LOAD_PROPERTY_SETTING, null, null);
                    
                } catch (FileNotFoundException ex) {
                    logger.error(ex.getMessage());
                }
            }
        });
        t1.start();

    }

    public void displayObservationSetting() {

        for (ObservationSet os : this.observationSetList) {
            System.out.println("observation \t\t\t" + os.getName() + "\n"
                    + "file \t\t\t" + os.getFile() + "\n"
                    + "species \t\t\t" + os.getSpecies() + "\n"
                    + "condition ionic strength [mM] \t" + os.getCondition().ionicStrength + "\n"
                    + "condition pH \t\t" + os.getCondition().pH + "\n"
                    + "model/prediction \t\t" + Boolean.toString(os.isActive()) + "/" + Boolean.toString(os.isPredict()) + "\n"
                    + "\n");
        }
    }

    public void loadDescriptorList(String s) {

        logger.info("Reading descriptor values from file " + s);
        System.out.println(ColorStatic.PURPLE + "Reading descriptor values from file " + s + ColorStatic.RESET);

        File sFile = new File(s);

        if (!sFile.exists()) {
            logger.warn("Descriptor value file " + sFile.getAbsolutePath() + " not found.");
        } else {

            ImportExport io = new ImportExport();
            this.projectDescriptor = io.importDescriptorList(sFile);

            changes.firePropertyChange(MoleculeStatic.CHANGE_CALC_DESCRIPTOR, null, null);
        }

    }

    public void loadModelSetting() {

    }

    /*
     calculates all descriptor entries as defined by descriptorSetList for all 
     observations as defined by structureSetList (if set active)
     */
    public void calcDescriptorAll() {

        if (descriptorSetList.isEmpty()) {
            System.out.println("No valid descriptor settings loaded.");
            return;
        }

        if (moleculeList.isEmpty()) {
            System.out.println("No valid structure information provided.");
            return;
        }

        Thread t2;
        t2 = new Thread(new Runnable() {

            @Override
            public void run() {

                if (settingProcessor == 1) {
                    long startTime = System.currentTimeMillis();
                    
                    try {
                        projectDescriptor.setDescriptorList(processDescriptorsParallel(moleculeList));

                        long stopTime = System.currentTimeMillis();
                        long elapsedTime = stopTime - startTime;
                        System.out.println(ColorStatic.BLUE + "Calculation time was " + elapsedTime / 1000 + " sec." + ColorStatic.RESET);
                        logger.info("Descriptors calculated in {} sec", elapsedTime / 1000);
                        changes.firePropertyChange(MoleculeStatic.CHANGE_CALC_DESCRIPTOR, null, null);

                    } catch (InterruptedException | ExecutionException ex) {
                        logger.error(ex.getMessage());
                    }

                } else {

                    processDescriptorsSingle(moleculeList);
                    exportDescriptorList();
                }

            }
        });
        t2.start();

    }

    public void exportDescriptorList() {

        File f = new File(this.projectPath + "/descriptorList.json");

        ImportExport io = new ImportExport();
        io.exportDescriptorListToJson(this.projectDescriptor, f);
    }

    public void processDescriptorsSingle(List<Molecule> inputs) {

        int countExist; 
        if(projectDescriptor.getDescriptorList() != null || !projectDescriptor.getDescriptorList().isEmpty()){
        countExist = projectDescriptor.getDescriptorList().size(); 
        System.out.println(ColorStatic.GREEN + "Resuming descriptor calculation for observation " + (countExist) + "." + ColorStatic.RESET) ;
        }else{
        countExist = 0; 
        System.out.println(ColorStatic.PURPLE + "No previously calculated descriptors loaded." + ColorStatic.RESET);
        }

        for (int i = countExist; i < inputs.size(); i++) {

            final Boolean active = observationSetList.get(i).isActive();

            // process your input here and compute the output
            List<DescriptorGroup> dg = new ArrayList<>();
            if (active == true) {
                try {
                    dg = calcDescriptorEntry(i);
                } catch (Exception e) {
                    System.out.println(ColorStatic.RED + e.getMessage() + ColorStatic.RESET);
                }
            }
            
            try {
                projectDescriptor.descriptorList.set(i, dg);
            } catch (Exception e) {
                projectDescriptor.descriptorList.add(dg);
            }
            changes.firePropertyChange(MoleculeStatic.CHANGE_CALC_DESCRIPTOR, null, null);
            exportDescriptorList();
        }

        changes.firePropertyChange(MoleculeStatic.CHANGE_CALC_DESCRIPTOR, null, null);

    }

    public List<List<DescriptorGroup>> processDescriptorsParallel(List<Molecule> inputs)
            throws InterruptedException, ExecutionException {

        int threads = Runtime.getRuntime().availableProcessors();
        ExecutorService service = Executors.newFixedThreadPool(threads);

        List<Future<List<DescriptorGroup>>> futures = new ArrayList<>();
        
        for (int i = 0; i < inputs.size(); i++) {

            final Boolean active = observationSetList.get(i).isActive();

            final int observationID = i;
            Callable<List<DescriptorGroup>> callable = new Callable<List<DescriptorGroup>>() {

                @Override
                public List<DescriptorGroup> call() throws Exception {

                    // process your input here and compute the output
                    if (active == true) {

                        return calcDescriptorEntry(observationID);
                    } else {
                        List<DescriptorGroup> emptyGroup = new ArrayList<>();
                        return emptyGroup;
                    }
                }
            };
            futures.add(service.submit(callable));
        }

        service.shutdown();

        List<List<DescriptorGroup>> outputs = new ArrayList<>();
        for (Future<List<DescriptorGroup>> future : futures) {
            outputs.add(future.get());
        }

        changes.firePropertyChange(MoleculeStatic.CHANGE_CALC_DESCRIPTOR, null, null);

        return outputs;
    }

    
    public List<DescriptorGroup> calcDescriptorEntry(int observationID) {

        List<DescriptorGroup> dGroup = new ArrayList<>();

        List<Double[]> sphereP = Vector.calcSphere(120);

        for (int i = 0; i < descriptorSetList.size(); i++) {
            logger.info("calculating desriptor set " + i);

            if (descriptorSetList.get(i).getProjection() != null) {

                
                if (descriptorSetList.get(i).getProjection().getSelectIO() == 1) {
                    sphereP = Vector.calcSphere(descriptorSetList.get(i).getProjection().orientation);

                    DescriptorGroup dg = this.calcDescriptorGroup(observationID, i, sphereP);
                    dGroup.add(dg);
                    sphereP = dg.getVector();

                } else {

                    System.out.println(ColorStatic.BLUE + "Using orientation identified by previous descriptor set." + ColorStatic.RESET);
                    DescriptorGroup dg = this.calcDescriptorGroup(observationID, i, sphereP);
                    dGroup.add(dg);
                }
            } else {
                logger.info("Descriptor set without orientation setting.");
                DescriptorGroup dg = this.calcDescriptorGroup(observationID, i, null);
                dGroup.add(dg);
            }
        }
        this.getMolecule(observationID).clearSurface();
        return dGroup;

    }

    
    public DescriptorGroup calcDescriptorGroup(int observationID, int descriptorID, List<Double[]> vector) {

        DescriptorSet descriptorSet = this.descriptorSetList.get(descriptorID);

        switch (descriptorSet.getName()) {
            case "DescrPlane":
                System.out.println(ColorStatic.BLUE + "Calculating plane descriptors." + ColorStatic.RESET);

                PlaneDescriptorGroup planeDescr = new PlaneDescriptorGroup(observationID, descriptorID);
                planeDescr.setVector(vector);
                try {
                    planeDescr.calcDescriptor();
                } catch (Exception e) {
                    System.out.println(e.getMessage());
                    System.out.println(ColorStatic.RED + "Could not calculate " + this.getDescriptorSet(descriptorID).getName()
                            + " for observation " + this.getObservationSet(observationID).getFile() + ColorStatic.RESET);
                    System.out.println(ColorStatic.RED + e.getMessage() + ColorStatic.RESET);
                }
                return planeDescr;

            case "DescrPatch":
                System.out.println(ColorStatic.BLUE + "Calculating patch descriptors." + ColorStatic.RESET);

                PatchDescriptorGroup patchDescr = new PatchDescriptorGroup(observationID, descriptorID);
                patchDescr.setVector(vector);
                try {
                    patchDescr.calcDescriptor();
                } catch (Exception e) {

                    System.out.println(ColorStatic.RED + "Could not calculate " + this.getDescriptorSet(descriptorID).getName()
                            + " for observation " + this.getObservationSet(observationID).getFile() + ColorStatic.RESET);
                    System.out.println(ColorStatic.RED + e.getMessage() + ColorStatic.RESET);
                }

                return patchDescr;

            case "DescrSphere":
                System.out.println(ColorStatic.BLUE + "Calculating sphere descriptors." + ColorStatic.RESET);

                SphereDescriptorGroup sphereDescr = new SphereDescriptorGroup(observationID, descriptorID);
                sphereDescr.setVector(vector);
                try {
                    sphereDescr.calcDescriptor();
                } catch (Exception e) {

                    System.out.println(ColorStatic.RED + "Could not calculate " + this.getDescriptorSet(descriptorID).getName()
                            + " for observation " + this.getObservationSet(observationID).getFile() + ColorStatic.RESET);
                    System.out.println(ColorStatic.RED + e.getMessage() + ColorStatic.RESET);
                }

                return sphereDescr;

            case "DescrSurface":
                System.out.println(ColorStatic.BLUE + "Calculating surface descriptors." + ColorStatic.RESET);

                SurfaceDescriptorGroup surfDescr = new SurfaceDescriptorGroup(observationID, descriptorID);
                try {
                    surfDescr.calcDescriptor();
                } catch (Exception e) {

                    System.out.println(ColorStatic.RED + "Could not calculate " + this.getDescriptorSet(descriptorID).getName()
                            + " for observation " + this.getObservationSet(observationID).getFile() + ColorStatic.RESET);
                    System.out.println(ColorStatic.RED + e.getMessage() + ColorStatic.RESET);
                }

                return surfDescr;

            case "DescrGeo":
                System.out.println(ColorStatic.BLUE + "Calculating shape descriptors." + ColorStatic.RESET);
                ShapeDescriptorGroup shapeDescr;
                shapeDescr = new ShapeDescriptorGroup(observationID, descriptorID);

                try {
                    shapeDescr.calcDescriptor();
                } catch (Exception e) {

                    System.out.println(ColorStatic.RED + "Could not calculate " + this.getDescriptorSet(descriptorID).getName()
                            + " for observation " + this.getObservationSet(observationID).getFile() + ColorStatic.RESET);
                    System.out.println(ColorStatic.RED + e.getMessage() + ColorStatic.RESET);
                }

                return shapeDescr;

            default:
                logger.error("Descriptor type could not be identified.");
                return null;
        }
    }

    public void showDescriptor() {
        int m = projectDescriptor.descriptorList.size();
        int n = projectDescriptor.descriptorList.get(0).size();
        int nD = 0;

        System.out.println(ColorStatic.RED + " Descriptors calculated: \n"
                + "observations: " + m + " \n "
                + "descriptor groups: " + n + " \n" + ColorStatic.RESET);

        for (int j = 0; j < n; j++) {
            System.out.print(j + ": ");

            List<Descriptor> dl = projectDescriptor.descriptorList.get(0).get(j).getDescriptor();
            for (Descriptor dl1 : dl) {
                nD++;
                System.out.print(" " + dl1.getValue());
            }
            System.out.println("");
        }
    }

    public int activeObservationCount() {
        int m = 0;
        for (ObservationSet observationSetList1 : observationSetList) {
            if (observationSetList1.isActive() == true) {
                m++;
            }
        }
        return m;
    }

    public int getObservationCount() {
        return this.observationSetList.size();
    }


    public int modelObservationCount() {
        int m = 0;
        for (ObservationSet observationSetList1 : observationSetList) {
            if (observationSetList1.isActive() == true & observationSetList1.isPredict() == false) {
                m++;
            }
        }
        return m;
    }

    public int predictObservationCount() {
        int m = 0;
        for (ObservationSet observationSetList1 : observationSetList) {
            if (observationSetList1.isActive() == true & observationSetList1.isPredict() == true) {
                m++;
            }
        }
        return m;
    }
    
    private void calcDescriptorMatrix() {

        if (projectDescriptor.descriptorList == null || projectDescriptor.descriptorList.isEmpty()) {
            this.descriptorMatrix = new Matrix(3, 3);
            return;
        }

        Thread t2 = new Thread(new Runnable() {

            @Override
            public void run() {

                int m = projectDescriptor.descriptorList.size();
                int nG = projectDescriptor.descriptorList.get(0).size();
                int nD = 0;

                for (int j = 0; j < nG; j++) {
                    DescriptorGroup dg = projectDescriptor.descriptorList.get(0).get(j);
                    nD = nD + dg.getDescriptor().size();
                }
                descriptorMatrix = new Matrix(m, nD);
                int ni;
                int mi = 0;

                for (int i = 0; i < m; i++) {
                    ObservationSet os = observationSetList.get(i);
                    List<DescriptorGroup> ds = projectDescriptor.descriptorList.get(i);
                    ni = 0;

                    for (DescriptorGroup dg : ds) {

                        for (Descriptor descriptor : dg.getDescriptor()) {
                            descriptorMatrix.data[mi][ni] = descriptor.getValue();
                            ni++;
                        }
                    }
                    mi++;
                }
            }
        });
        t2.start();
        try {
            t2.join();
        } catch (InterruptedException ex) {
            logger.error(ex.getMessage());
        }
    }

    private void calcModelDescriptorMatrix() {

        if (projectDescriptor.descriptorList == null || projectDescriptor.descriptorList.isEmpty()) {
            System.out.println(ColorStatic.RED + "Error: projectDescriptor.descriptorList is Empty.");
            return;
        }

        Thread t2 = new Thread(new Runnable() {

            @Override
            public void run() {

                int m = modelObservationCount();
                int nG = projectDescriptor.descriptorList.get(0).size();
                int nD = 0;

                for (int j = 0; j < nG; j++) {
                    DescriptorGroup dg = projectDescriptor.descriptorList.get(0).get(j);
                    nD = nD + dg.getDescriptor().size();
                }

                modelDescriptorMatrix = new Matrix(m, nD);
                int ni;
                int mi = 0;

                for (int i = 0; i < observationSetList.size(); i++) {
                    ObservationSet os = observationSetList.get(i);
                    List<DescriptorGroup> ds = projectDescriptor.descriptorList.get(i);

                    if (os.isActive() == true & os.isPredict() == false) {

                        ni = 0;

                        for (DescriptorGroup dg : ds) {

                            for (Descriptor descriptor : dg.getDescriptor()) {
                                modelDescriptorMatrix.data[mi][ni] = descriptor.getValue();
                                ni++;
                            }
                        }
                        mi++;
                    }
                }
            }
        });
        t2.start();
        try {
            t2.join();
        } catch (InterruptedException ex) {
            logger.error(ex.getMessage());
        }
    }

    public Matrix getModelDescriptorMatrix() {

        this.calcModelDescriptorMatrix();
        return modelDescriptorMatrix;

    }

    public Matrix getDescriptorMatrix() {
        this.calcDescriptorMatrix();
        return this.descriptorMatrix;
    }

    public Matrix getPredictDescriptorMatrix() {
        this.calcPredictDescriptorMatrix();

        return predictDescriptorMatrix;
    }

    public DescriptorGroup getDescriptorGroup(int molPos, int descrPos) {

        DescriptorGroup dG = null;
        List<DescriptorGroup> descrGroupList = projectDescriptor.descriptorList.get(molPos);

        int pos = descrPos;

        for (DescriptorGroup descrGroupList1 : descrGroupList) {

            if (pos < descrGroupList1.descriptorList.size()) {
                dG = descrGroupList1;

                break;
            } else {
                pos = pos - descrGroupList1.descriptorList.size();
            }
        }
        return dG;
    }

    public int getDescriptorPositionInGroup(int molPos, int descrPos) {

        List<DescriptorGroup> descrGroupList = projectDescriptor.descriptorList.get(molPos);

        int pos = descrPos;

        for (DescriptorGroup descrGroupList1 : descrGroupList) {

            if (pos < descrGroupList1.descriptorList.size()) {
                break;
            } else {
                pos = pos - descrGroupList1.descriptorList.size();
            }
        }
        return pos;
    }

    private void calcPredictDescriptorMatrix() {

        Thread t2 = new Thread(new Runnable() {

            @Override
            public void run() {

                int m = predictObservationCount();
                int nG = projectDescriptor.descriptorList.get(0).size();
                int nD = 0;

                for (int j = 0; j < nG; j++) {
                    DescriptorGroup dg = projectDescriptor.descriptorList.get(0).get(j);
                    nD = nD + dg.getDescriptor().size();
                }

                predictDescriptorMatrix = new Matrix(m, nD);

                int ni;
                int mi = 0;
                for (int i = 0; i < observationSetList.size(); i++) {
                    ObservationSet os = observationSetList.get(i);
                    List<DescriptorGroup> ds = projectDescriptor.descriptorList.get(i);

                    if (os.isActive() == true & os.isPredict() == true) {

                        ni = 0;
                        for (DescriptorGroup dg : ds) {

                            for (Descriptor descriptor : dg.descriptorList) {
                                predictDescriptorMatrix.data[mi][ni] = descriptor.getValue();
                                ni++;
                            }
                        }
                        mi++;
                    }
                }
            }
        });

        t2.start();

        try {
            t2.join();
        } catch (InterruptedException ex) {
            logger.error(ex.getMessage());
        }
    }

    public List<ObservationSet> getObservationSetList() {
        return this.observationSetList;
    }

    public List<Double> getModelObservationValues() {
        List<Double> ob = new ArrayList<>();

        for (ObservationSet os : this.observationSetList) {
            if (os.isActive() == true & os.isPredict() == false) {
                ob.add(os.getResponse());
            }
        }
        return ob;
    }

    public List<ObservationSet> getModelObservationSetList() {
        List<ObservationSet> ob = new ArrayList<>();

        for (ObservationSet os : this.observationSetList) {
            if (os.isActive() == true & os.isPredict() == false) {
                ob.add(os);
            }
        }
        return ob;
    }

    public List<Double> getPredictObservationValues() {
        List<Double> ob = new ArrayList<>();

        for (ObservationSet os : this.observationSetList) {
            if (os.isActive() == true & os.isPredict() == true) {
                ob.add(os.getResponse());
            }
        }
        return ob;
    }

    public List<ObservationSet> getPredictObservationSetList() {
        List<ObservationSet> ob = new ArrayList<>();

        for (ObservationSet os : this.observationSetList) {
            if (os.isActive() == true & os.isPredict() == true) {
                ob.add(os);
            }
        }
        return ob;
    }

    public Molecule getMolecule(int i) {
        return moleculeList.get(i);
    }

    public DescriptorSet getDescriptorSet(int i) {
        return this.descriptorSetList.get(i);
    }

    public ObservationSet getObservationSet(int i) {
        return this.observationSetList.get(i);
    }

    public int getMoleculeListSize() {
        return moleculeList.size();
    }

    public void setMolecule(int i, Molecule structureSet) {
        try {
            this.moleculeList.set(i, structureSet);
        } catch (IndexOutOfBoundsException ex) {
            logger.error(ex.getMessage());
        }
    }

    public List<String> getDescriptorName() {

        List<String> descriptorName = new ArrayList<>();
        //
        if (projectDescriptor.descriptorList == null || projectDescriptor.descriptorList.isEmpty()) {
            return null;
        }

        List<DescriptorGroup> descriptorGroupList = projectDescriptor.descriptorList.get(0);
        descriptorGroupList.stream().forEach((dg) -> {
            dg.descriptorList.stream().forEach((descriptor) -> {
                descriptorName.add(descriptor.getName());
            });
        });
        return descriptorName;
    }

    public String getProjectPath() {
        return projectPath;
    }

    public void setProjectPath(String projectPath) {
        this.projectPath = projectPath;
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        changes.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
        changes.removePropertyChangeListener(listener);
    }
}
