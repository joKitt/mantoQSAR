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
package org.mantoQSAR.core.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import org.mantoQSAR.core.math.MathUtil;
import org.mantoQSAR.core.math.Matrix;
import org.mantoQSAR.core.util.ColorStatic;

public class ReplacementModel extends ModelClass {

    private final double rcondUpper = 1e-25;

    List<double[]> vecTOT;
    List<Double> errTOT;

    public int[] vecSingle = null;
    double errSingle;
    Boolean[] descrExclude;
    List<List<double[]>> TOT;

    int numDesc;

    double s_err;

    public ReplacementModel() {

        super();
        initReplacementMethod();

    }

    private void initReplacementMethod() {

        // set for development purpose
        numDesc = 7;

        this.vecTOT = new ArrayList<>();
        this.errTOT = new ArrayList<>();
        this.TOT = new ArrayList<>();

        this.errSingle = 1000;
    }

    @Override
    public void calcModel() {

        getEmptyDescriptor();

        Random rng = new Random();
        int r_m = getModelDescriptor().getColumnCount() - 1;
        Set<Integer> vecSet = new LinkedHashSet<>();

        while (vecSet.size() < numDesc) {
            Integer next = rng.nextInt(r_m);

            if (descrExclude[next] == false) {
                vecSet.add(next);
            }
        }

        Integer[] VecI = new Integer[vecSet.size()];
        int[] vecI = new int[vecSet.size()];

        VecI = vecSet.toArray(VecI);

        for (int i = 0; i < VecI.length; i++) {
            vecI[i] = VecI[i];
        }
        erm(vecI);

        System.out.println("model calculated");
    }

    /*
     returns a List of boolean indicating descriptors to be considered in model 
     generation; descriptors taking the value zero for all observations are 
     excluded (true); 
     */
    private void getEmptyDescriptor() {

        Matrix modelDescr = this.getModelDescriptor();
        this.descrExclude = new Boolean[modelDescr.getColumnCount()];

        for (int i = 0; i < modelDescr.getColumnCount(); i++) {

            this.descrExclude[i] = true;

            for (int j = 0; j < modelDescr.getRowCount(); j++) {

                double a = modelDescr.data[j][i];

                if (a > 0.00000001 || a < -0.00000001) {
                    this.descrExclude[i] = false;
                }
            }
        }
    }

    /* Inverse Enhanced Replacement Method (IERM) algorithm according to 
     Mercader et al. 
    
     */
    private void ierm(int numDesc, int[] vec) {

        Integer r_m = this.getModelDescriptor().getColumnCount();
        getEmptyDescriptor();

        if (r_m < numDesc) {
            logger.error("Not enough descriptors calculated. \n" + r_m + " descriptors counted; model requires " + numDesc);
            return;
        }

        logger.info("start ierm calculation with " + numDesc + " descriptors");

        int[] vecI = vec;
        int ii = 0;
        int c = 0;
        while (c < numDesc) {
            if (this.descrExclude[ii] == false) {
                vecI[c] = ii;
                c++;
            }
            ii++;
        }

        double lindep = ld(vecI);
        Random rng = new Random();

        if (lindep == 100) {
            for (int j = 0; j < 2e31 - 1; j++) {

                Set<Integer> vecSet = new LinkedHashSet<>();

                while (vecSet.size() < numDesc) {
                    Integer next = rng.nextInt(r_m);

                    if (this.descrExclude[next] == false) {
                        vecSet.add(next);
                    }
                }

                Integer[] VecI = new Integer[vecSet.size()];
                VecI = vecSet.toArray(VecI);

                for (int i = 0; i < VecI.length; i++) {
                    vecI[i] = VecI[i];
                }
                lindep = ld(vecI);
                if (lindep == 0.0) {
                    break;
                }
            }

        }
        this.rmt_inv(vecI);

        erm(this.vecSingle);

    }

    /* Enhanced Replacement Method (ERM) algorithm according to 
     Mercader et al.  */
    private void erm(int[] vec) {

        logger.info("start vector rmt: " + Arrays.toString(vec));

        getEmptyDescriptor();

        List<List<double[]>> sTOT = new ArrayList<>();

        vecTOT.clear();
        errTOT.clear();

        int pos;
        int n_v = vec.length;

        double[] coeff;

        List<Integer> Po = new ArrayList<>();

        for (int i = 0; i < n_v; i++) {

            double sr = this.rms(vec);
            double[] a = new double[vec.length + 1];
            a[0] = sr;

            for (int j = 0; j < vec.length; j++) {
                a[j + 1] = vec[j];
            }

            List<double[]> A = new ArrayList<>();
            A.add(a);
            sTOT.add(A);

            double[] vecA = rmsr(vec, i);

            if (!Po.isEmpty()) {
                Po.set(0, i);
            } else {
                Po.add(i);
            }

            if (n_v == 1) {
                logger.info("target vector is 1.");
                double sVecTot[] = new double[vecA.length + 1];
                sVecTot[0] = 1;
                System.arraycopy(vecA, 0, sVecTot, 1, vecA.length);

                sTOT.clear();
                A = new ArrayList<>();
                A.add(vecA);
                sTOT.add(A);
                return;
            }

            int[] vecI = new int[vecA.length - 1];
            for (int k = 0; k < vecA.length - 1; k++) {
                vecI[k] = (int) vecA[k + 1];
            }

            coeff = rmder(vecI, false);

            pos = getPositionMaxCoeff(coeff, Po);
            Po.add(pos);

            logger.trace("Po value " + Po.toString());

            sTOT.get(i).add(vecA);

            for (int j = 2; j < n_v; j++) {

                logger.trace("iteration loop 1: " + j);
                vecA = rmsr(vecI, pos);
                vecI = new int[vecA.length - 1];
                for (int k = 0; k < vecA.length - 1; k++) {
                    vecI[k] = (int) vecA[k + 1];
                }

                sTOT.get(i).add(vecA);

                if (j + 1 == n_v) {
                    Po.clear();
                    break;
                }

                coeff = rmder(vecI);
                pos = getPositionMaxCoeff(coeff, Po);
                Po.add(pos);

                logger.trace("testing target vector " + Po.toString());
            }

            int j;

            for (j = 1; j < 100; j++) {

                logger.trace("iteration loop 2: " + j);
                pos = getPositionMaxCoeff(coeff, Po);

                if (!Po.isEmpty()) {
                    Po.set(0, pos);
                } else {
                    Po.add(pos);
                }

                for (int k = 1; k <= n_v; k++) {

                    logger.trace("loop A2_A iteration " + k);

                    vecA = rmsr(vecI, pos);
                    sTOT.get(i).add(vecA);

                    vecI = new int[vecA.length - 1];
                    for (int m = 0; m < vecA.length - 1; m++) {
                        vecI[m] = (int) vecA[m + 1];
                    }

                    if (k == n_v) {
                        Po.clear();
                        break;
                    }

                    coeff = rmder(vecI);
                    pos = getPositionMaxCoeff(coeff, Po);
                    Po.add(pos);

                    logger.trace("testing target vector " + Po.toString());
                }

                List<double[]> subList = sTOT.get(i);

                if ((j > 3) & (subList.size() > ((2 * n_v) + 1))) {
                    double[] listLast = subList.get(subList.size() - 1);
                    double[] listComp = subList.get(subList.size() - (2 * n_v));

                    if (Arrays.equals(listLast, listComp)) {
                        logger.info("erm primary algorithm converged");
                        break;
                    }
                }
            }

            int jj;
            for (jj = j; jj < j + 100; jj++) {

                logger.trace("iteration loop 2: " + j);
                pos = getPositionMaxCoeff(coeff, Po);

                if (!Po.isEmpty()) {
                    Po.set(0, pos);
                } else {
                    Po.add(pos);
                }

                for (int k = 0; k < n_v; k++) {

                    vecA = rma2(vecI, pos);
                    sTOT.get(i).add(vecA);
                    vecI = new int[vecA.length - 1];
                    for (int m = 0; m < vecA.length - 1; m++) {
                        vecI[m] = (int) vecA[m + 1];
                    }

                    if (k + 1 == n_v) {
                        Po.clear();
                    }
                    coeff = rmder(vecI);
                    pos = getPositionMaxCoeff(coeff, Po);
                    Po.add(pos);
                }

                List<double[]> subList = sTOT.get(i);

                if ((subList.size() > ((4 * n_v) + 1))) {
                    double[] listLast = subList.get(subList.size() - 1);
                    double[] listComp = subList.get(subList.size() - (4 * n_v));

                    if (Arrays.equals(listLast, listComp)) {
                        logger.info("erm secondary algorithm converged");
                        break;
                    }
                }
            }

            List<double[]> subList = sTOT.get(i);

            double val = 1000;
            double[] vecQ = new double[n_v + 1];
            for (int k = 0; k < subList.size(); k++) {

                if (val > subList.get(k)[0]) {
                    val = subList.get(i)[0];
                    vecQ = subList.get(i);
                }
            }
            vecI = new int[vecQ.length - 1];
            for (int m = 0; m < vecQ.length - 1; m++) {
                vecI[m] = (int) vecQ[m + 1];
            }

            for (int jjj = jj; jjj < jj + 100; jjj++) {

                logger.trace("iteration loop 2: " + j);
                pos = getPositionMaxCoeff(coeff, Po);

                if (!Po.isEmpty()) {
                    Po.set(0, pos);
                } else {
                    Po.add(pos);
                }

                for (int k = 0; k < n_v; k++) {

                    vecA = rmsr(vecI, pos);
                    sTOT.get(i).add(vecA);
                    vecI = new int[vecA.length - 1];
                    for (int m = 0; m < vecA.length - 1; m++) {
                        vecI[m] = (int) vecA[m + 1];
                    }

                    if (k + 1 == n_v) {
                        Po.clear();
                    }
                    coeff = rmder(vecI);
                    pos = getPositionMaxCoeff(coeff, Po);
                    Po.add(pos);
                }

                subList = sTOT.get(i);

                if ((subList.size() > ((2 * n_v) + 1))) {
                    double[] listLast = subList.get(subList.size() - 1);
                    double[] listComp = subList.get(subList.size() - (2 * n_v));

                    if (Arrays.equals(listLast, listComp)) {
                        logger.info("erm tertiary algorithm converged");
                        break;
                    }
                }
            }

            double vecErr = 1000.0;
            double[] vecP = new double[vec.length + 1];
            subList = sTOT.get(i);

            for (double[] subList1 : subList) {
                if (subList1[0] < vecErr) {
                    vecErr = subList1[0];
                    System.arraycopy(subList1, 0, vecP, 0, subList1.length);
                    logger.trace("Update vecErr " + vecErr);
                }
            }

            if (vecErr < this.errSingle) {
                this.errSingle = vecErr;
                this.vecSingle = new int[n_v];
                for (int k = 0; k < n_v; k++) {
                    this.vecSingle[k] = (int) vecP[k + 1];
                }
            }

            vecTOT.add(vecP);
            errTOT.add(vecErr);

        }
        logger.info("end vector " + Arrays.toString(vecSingle));

    }

    private int getPositionMinCoeff(double[] coeff, List<Integer> Po) {

        if (coeff == null) {
            logger.error("coeff is null in getPositionMinCoeff");
            return 0;
        }

        if (Po == null) {
            throw new NullPointerException("Po is null");
        }

        List<Double> Coer = new ArrayList<>();
        boolean excl;

        for (int k = 0; k < coeff.length; k++) {

            excl = false;
            for (Integer Po1 : Po) {
                if (k == (int) Po1) {
                    excl = true;
                }
            }
            if (excl == false) {
                Coer.add(coeff[k]);
            }
        }

        double minVal = Double.MAX_VALUE;
        for (Double Coer1 : Coer) {
            if (Coer1 < minVal) {
                minVal = Coer1;
            }
        }

        int pos;
        for (int k = 0; k < coeff.length; k++) {
            if (coeff[k] == minVal) {
                pos = k;
                return pos;
            }
        }

        logger.error("Illegal value exception");
        logger.error("Coeff " + Arrays.toString(coeff));
        logger.error("PO: " + Po.toString());

        return 0;
    }

    private int getPositionMaxCoeff(double[] coeff, List<Integer> Po) {

        List<Double> Coer = new ArrayList<>();

        boolean excl;
        for (int k = 0; k < coeff.length; k++) {

            excl = false;
            for (Integer Po1 : Po) {
                if (k == (int) Po1) {
                    excl = true;
                }
            }
            if (excl == false) {
                Coer.add(coeff[k]);
            }
        }

        double maxVal = Double.MIN_VALUE;
        for (Double Coer1 : Coer) {
            if (Coer1 > maxVal) {
                maxVal = Coer1;
            }
        }

        int pos;
        for (int k = 0; k < coeff.length; k++) {
            if (coeff[k] == maxVal) {
                pos = k;
                return pos;
            }
        }

        logger.error("Illegal value exception");
        logger.error("Coeff " + Arrays.toString(coeff));
        logger.error("PO: " + Po.toString());
        return 0;

    }

    private void rmt(int[] vec) {
        this.rmt(vec, false);
    }

    private void rmt_inv(int[] vec) {
        this.rmt(vec, true);
    }

    private void rmt(int[] vec, Boolean inv) {

        Boolean inverse = false;

        if (inv != null) {
            inverse = inv;
        }

        logger.info("start rmt calculation with ");
        logger.info("start vector: " + Arrays.toString(vec));

        List<List<double[]>> sTOT = new ArrayList<>();
        List<double[]> sVecTOT = new ArrayList<>();

        vecTOT.clear();
        errTOT.clear();

        int n_v = vec.length;
        double[] vecA = new double[n_v + 1];
        double[] coeff;
        double[] coeffTemp;
        int pos;

        List<Integer> Po = new ArrayList<>();

        for (int i = 0; i < n_v; i++) {

            double sr;
            if (inverse == false) {
                sr = this.rms_inv(vec);
            } else {
                sr = this.rms(vec);
            }

            double[] a = new double[vec.length + 1];
            a[0] = sr;

            for (int j = 0; j < vec.length; j++) {
                a[j + 1] = vec[j];
            }

            List<double[]> A = new ArrayList<>();
            A.add(a);
            sTOT.add(A);

            if (inverse == false) {
                vecA = this.rmsr(vec, i);
            } else {
                vecA = this.rmsr_inv(vec, i);
            }

            if (!Po.isEmpty()) {
                Po.set(0, i);
            } else {
                Po.add(i);
            }

            if (n_v == 1) {

                double sVecTot[] = new double[vecA.length + 1];
                sVecTot[0] = 1;
                System.arraycopy(vecA, 0, sVecTot, 1, vecA.length);

                sTOT.clear();
                A = new ArrayList<>();
                A.add(vecA);
                sTOT.add(A);
                return;
            }

            int[] vecI = new int[vecA.length - 1];
            for (int k = 0; k < vecA.length - 1; k++) {
                vecI[k] = (int) vecA[k + 1];
            }

            coeff = this.rmder(vecI, false);

            if (inverse == false) {
                pos = getPositionMaxCoeff(coeff, Po);
            } else {
                pos = getPositionMinCoeff(coeff, Po);
            }
            Po.add(pos);

            sTOT.get(i).add(vecA);

            for (int j = 2; j < n_v; j++) {

                logger.trace("iteration loop 1: " + j);

                if (inverse == false) {
                    vecA = this.rmsr(vec, i);
                } else {
                    vecA = this.rmsr_inv(vec, i);
                }

                vecI = new int[vecA.length - 1];
                for (int k = 0; k < vecA.length - 1; k++) {
                    vecI[k] = (int) vecA[k + 1];
                }

                sTOT.get(i).add(vecA);

                if (j + 1 == n_v) {
                    Po.clear();
                    break;
                }

                coeffTemp = this.rmder(vecI, inverse);
                if (coeffTemp != null) {
                    coeff = coeffTemp;
                    if (inverse == false) {
                        pos = getPositionMaxCoeff(coeff, Po);
                    } else {
                        pos = getPositionMinCoeff(coeff, Po);
                    }
                    Po.add(pos);
                } else {
                    logger.error("coeff was null");
                }
                logger.trace("testing target vector " + Po.toString());
            }

            for (int j = 1; j < 100; j++) {

                if (inverse == false) {
                    pos = getPositionMaxCoeff(coeff, Po);
                } else {
                    pos = getPositionMinCoeff(coeff, Po);
                }

                if (!Po.isEmpty()) {
                    Po.set(0, pos);
                } else {
                    Po.add(pos);
                }

                for (int k = 1; k <= n_v; k++) {

                    if (inverse == false) {
                        vecA = this.rmsr(vec, i);
                    } else {
                        vecA = this.rmsr_inv(vec, i);
                    }

                    sTOT.get(i).add(vecA);
                    vecI = new int[vecA.length - 1];
                    for (int m = 0; m < vecA.length - 1; m++) {
                        vecI[m] = (int) vecA[m + 1];
                    }

                    if (k == n_v) {
                        Po.clear();
                        break;
                    }

                    coeffTemp = this.rmder(vecI, inverse);

                    if (coeffTemp != null) {
                        coeff = coeffTemp;
                    }
                    if (inverse == false) {
                        pos = this.getPositionMaxCoeff(coeff, Po);
                    } else {
                        pos = this.getPositionMinCoeff(coeff, Po);
                    }
                    Po.add(pos);
                }

                List<double[]> subList = sTOT.get(i);

                if ((j > 3) & (subList.size() > ((2 * n_v) + 1))) {
                    double[] listLast = subList.get(subList.size() - 1);
                    double[] listComp = subList.get(subList.size() - (2 * n_v));

                    if (Arrays.equals(listLast, listComp)) {
                        Po.clear();
                        logger.info("algorithm converged");
                        break;
                    }
                }
            }

            double vecErr;
            if (inverse) {
                vecErr = Double.MIN_VALUE;
            } else {
                vecErr = Double.MAX_VALUE;
            }

            double[] vecP = new double[vec.length + 1];
            List<double[]> subList = sTOT.get(i);

            for (double[] subList1 : subList) {
                if (inverse) {
                    if (subList1[0] > vecErr) {
                        vecErr = subList1[0];
                        System.arraycopy(subList1, 0, vecP, 0, subList1.length);
                        this.vecSingle = new int[n_v];
                        for (int k = 0; k < n_v; k++) {
                            this.vecSingle[k] = (int) vecP[k + 1];
                        }
                    }
                } else {
                    if (subList1[0] < vecErr) {
                        vecErr = subList1[0];
                        System.arraycopy(subList1, 0, vecP, 0, subList1.length);
                        this.vecSingle = new int[n_v];
                        for (int k = 0; k < n_v; k++) {
                            this.vecSingle[k] = (int) vecP[k + 1];
                        }
                    }
                }
            }
            if (vecErr < errSingle) {
                errSingle = vecErr;
                this.vecSingle = new int[n_v];
                for (int k = 0; k < n_v; k++) {
                    this.vecSingle[k] = (int) vecP[k + 1];
                }
            }
            vecTOT.add(vecP);
            errTOT.add(vecErr);
        }
    }

    public void calcStepwise(final int maxDesc) {

        Thread t1 = new Thread(new Runnable() {

            @Override
            public void run() {
                List<double[]> stepResult = stepwise(maxDesc);

                System.out.println(ColorStatic.GREEN + "Stepwise results" + ColorStatic.RESET);
                for (double[] stepResult1 : stepResult) {
                    System.out.println("\t \t" + Arrays.toString(stepResult1));
                }
            }
        });
        t1.start();
    }

    private Matrix calcCoefficients() {

        if (vecSingle == null) {
            return null;
        }

        int[] vec = this.vecSingle;
        logger.debug("vec single " + Arrays.toString(this.vecSingle));

        List<Double> P = this.getModelProperty();
        Matrix X = this.getModelDescriptor().getColumns(vec);

        int m = X.data.length;
        int n = (X.getColumnCount() + 1);

        Matrix XX = new Matrix(m, n);
        for (int i = 0; i < m; i++) {
            XX.data[i][0] = 1.0;
            for (int j = 0; j < X.getColumnCount(); j++) {
                XX.data[i][j + 1] = X.getData(i, j);
            }
        }

        Matrix XXt = XX.transpose();
        Matrix XTXI = XXt.times(XX);
        Matrix XTXIinv = XTXI.invert();

        Matrix PP = new Matrix(P.size(), 1);
        for (int i = 0; i < P.size(); i++) {
            PP.data[i][0] = P.get(i);
        }

        Matrix a = XTXIinv.times(XXt);
        Matrix coef = a.times(PP);

        return coef;
    }

    /* returns List of double values representing prediction values from model 
     for descriptors in Matrix descriptorMatrix
     */
    public List<Double> getModelResponse(Matrix descriptorMatrix) {

        if (this.vecSingle == null) {
            return null;
        }

        Matrix xP = descriptorMatrix.getColumns(this.vecSingle);
        int m = xP.getRowCount();
        int n = (xP.getColumnCount() + 1);

        Matrix XX = new Matrix(m, n);
        for (int i = 0; i < m; i++) {

            XX.data[i][0] = 1.0;
            for (int j = 0; j < xP.getColumnCount(); j++) {
                XX.data[i][j + 1] = xP.getData(i, j);
            }
        }

        Matrix coef = this.calcCoefficients();
        List<Double> response = new ArrayList<>();

        for (int i = 0; i < m; i++) {
            Double d = 0.0;
            for (int j = 0; j < n; j++) {
                Double a = XX.data[i][j] * coef.data[j][0];
                XX.data[i][j] = a;
                d = d + a;
            }
            response.add(d);
        }
        logger.info("calculated prediction: " + response.toString());
        return response;
    }

    public List<Double> getModelResponse() {
        return this.getModelResponse(this.getModelDescriptor());
    }

    public List<Integer> getDescriptorsSelected() {
        if (vecSingle == null) {
            return null;
        }
        List<Integer> descrList = new ArrayList<>(this.vecSingle.length);
        for (int i = 0; i < this.vecSingle.length; i++) {
            descrList.add(this.vecSingle[i]);
        }
        return descrList;
    }

    public List<Double> getCoefficients() {

        if (this.vecSingle == null) {
            return null;
        }
        Matrix coef = this.calcCoefficients();
        List<Double> coefList = new ArrayList<>(coef.getRowCount());

        for (int i = 0; i < coef.data.length; i++) {
            coefList.add(coef.data[i][0]);
        }
        return coefList;
    }

    public List<Double> getPredResponse() {

        return this.getModelResponse(this.getPredictDescriptor());
    }

    private List<double[]> stepwise(int maxDesc) {
        List<double[]> resSW = new ArrayList<>();

        int n_m = this.getModelDescriptor().getColumnCount();
        double smin;

        List<Integer> VecJ = new ArrayList<>();
        int[] vecK = null;

        logger.debug("stepwise iteration to " + maxDesc + " descriptor model");

        for (int i = 1; i <= maxDesc; i++) {
            smin = Double.MAX_VALUE;

            for (int j = 0; j < n_m; j++) {
                int[] vecI = new int[VecJ.size() + 1];
                for (int k = 0; k < VecJ.size(); k++) {
                    vecI[k] = VecJ.get(k);
                }
                vecI[vecI.length - 1] = j;

                double serr = rms(vecI);

                if (serr < smin) {
                    smin = serr;
                    vecK = new int[vecI.length];
                    System.arraycopy(vecI, 0, vecK, 0, vecI.length);
                }
            }

            VecJ.clear();
            for (int k = 0; k < vecK.length; k++) {
                VecJ.add(vecK[k]);
            }

            double[] A = new double[maxDesc + 1];
            A[0] = smin;
            for (int k = 0; k < maxDesc; k++) {
                try {
                    A[k + 1] = vecK[k];
                } catch (Exception e) {
                    A[k + 1] = Double.NaN;
                }
            }
            logger.debug(Arrays.toString(A));
            System.out.println(Arrays.toString(A));
            resSW.add(A);
        }
        return resSW;
    }

    /* ld tests the lineal dependence of the descriptors in vec including 
     the regression constant
     */
    private Double ld(int[] vec) {

        Double lindep = 0.0;

        int n = this.getModelProperty().size();
        int k = vec.length;
        Matrix X = this.getModelDescriptor().getColumns(vec);

        double[][] xx = new double[X.getRowCount()][X.getColumnCount() + 1];

        for (int i = 0; i < X.getRowCount(); i++) {

            xx[i][0] = 1.0;
            for (int j = 0; j < X.getColumnCount(); j++) {
                xx[i][j + 1] = X.getData(i, j);
            }
        }

        Matrix XX = new Matrix(xx);

        int nu = k + 1;
        if (XX.rank() < nu) {
            lindep = 100.0;
        }
        return lindep;
    }

    /* replaces all descriptors from Mat except the ones in the given descriptor 
     vector in the given position Pth and returns the vector with lower S_err,
     using multiple linear regression analysis. 
     (The initial descriptor is changed even if it has the lower S_err) */
    private double[] rma2(int[] vec, Integer pth) {

        Matrix modelDescr = this.getModelDescriptor();
        int k = modelDescr.getRowCount();
        int n_m = modelDescr.getColumnCount();

        List<Integer> num = new ArrayList<>();

        List<Integer> Vec = new ArrayList<>();
        for (int i = 0; i < vec.length; i++) {
            Vec.add(vec[i]);
        }

        for (int i = 0; i < n_m; i++) {

            if (Vec.contains(i)) {
            } else {
                num.add(i);
            }
        }

        int n_n = num.size();
        int desc = 1;

        double smin = 10000000;
        for (int j = 0; j < n_n; j++) {
            vec[pth] = num.get(j);

            double serr = rms(vec);
            if (serr < smin) {
                smin = serr;
                desc = num.get(j);
            }
        }

        vec[pth] = desc;
        double[] a = new double[vec.length + 1];
        a[0] = smin;

        for (int j = 0; j < vec.length; j++) {
            a[j + 1] = vec[j];
        }
        return a;
    }

    /* rmder returns a vector containing the relative standard deviation of the
     regression coefficients. The constant of the regression is not included. */
    private double[] rmder(int[] vec) {
        return rmder(vec, false);
    }

    private double[] rmder(int[] vec, Boolean inv) {

        Boolean inverse = false;

        if (inv != null) {
            inverse = inv;
        }
        int n = this.getModelProperty().size();
        int k = vec.length;
        double[] coeff = new double[k];

        // prime with high number
        for (int i = 0; i < k; i++) {
            coeff[i] = 10000.0;
        }

        Matrix X = this.getModelDescriptor().getColumns(vec);

        // Check that independent (X) and dependent (P) data have compatible dimensions
        if (n != X.getRowCount()) {
            throw new RuntimeException("rmder() - Illegal vector length.");
        }

        // Solve for the regression coefficients using ordinary least-squares
        double[][] xx = new double[X.getRowCount()][X.getColumnCount() + 1];

        for (int i = 0; i < X.getRowCount(); i++) {
            xx[i][0] = 1.0;
            for (int j = 0; j < X.getColumnCount(); j++) {
                xx[i][j + 1] = X.getData(i, j);
            }
        }

        Matrix XX = new Matrix(xx);
        Matrix XXt = XX.transpose();
        Matrix XTXI = XXt.times(XX);
        Matrix XTXIinv = XTXI.invert();

        // get reciprocal condition of matrix 
        double rcond = 1 / (XTXIinv.norm1() * XTXI.norm1());

        if (inverse == true) {
            if (rcond < this.rcondUpper || Double.isNaN(rcond)) {
                this.s_err = 1e-26;
                logger.debug("return for rcond out of bounds");
                return coeff;
            }
        } else {
            if (rcond < this.rcondUpper || Double.isNaN(rcond)) {
                this.s_err = 1000000;
                logger.debug("return for rcond out of bounds");
                return coeff;
            }
        }

        List<Double> modelProp = this.getModelProperty();
        int length = modelProp.size();
        double[][] p = new double[length][1];
        for (int i = 0; i < length; i++) {
            p[i][0] = modelProp.get(i);
        }

        Matrix PP = new Matrix(p);
        Matrix a = XTXIinv.times(XXt);
        Matrix Coefficients = a.times(PP);

        for (int i = 0; i < Coefficients.getRowCount(); i++) {
            for (int j = 0; j < Coefficients.getColumnCount(); j++) {

                if (Double.isNaN(Coefficients.data[i][j])) {
                    Coefficients.data[i][j] = 1.0;
                }
            }
        }

        Matrix PP_hat = XX.times(Coefficients);
        Matrix residuals = PP.minus(PP_hat);
        Matrix b = residuals.transpose().times(residuals);
        this.s_err = Math.sqrt(b.getData(0, 0) / (n - k - 1));

        if (Double.isInfinite(this.s_err) || Double.isNaN(this.s_err)) {
            if (inverse == true) {
                this.s_err = 1e-26;
            } else {
                this.s_err = 1000000;
            }
        }

        Matrix covariance = XTXIinv.times(Math.pow(s_err, 2));

        double[] C = covariance.getDiagonal();
        double[] Er = new double[C.length];
        double c;

        for (int i = 0; i < C.length; i++) {
            c = Math.sqrt(C[i]);
            Er[i] = Math.abs(c / Coefficients.data[i][0]) * 100;
        }

        for (int i = 0; i < vec.length; i++) {
            if (Double.isNaN(Er[i + 1])) {
                coeff[i] = 9.9;
            } else {
                coeff[i] = Er[i + 1];
            }
        }
        return coeff;
    }

    /* rms evaluates selection of descriptors using multiple linear regression
     analysis.  */
    private Double rms(int[] vec) {
        return rms(vec, false);
    }

    private Double rms(int[] vec, Boolean inv) {

        Boolean inverse = false;

        List<Double> modelProp = this.getModelProperty();
        Matrix X = this.getModelDescriptor().getColumns(vec);

        int lengthP = modelProp.size();

        if (inv != null) {
            inverse = inv;
        }

        int n = modelProp.size();
        int k = vec.length;

        // Check that independent (X) and dependent (P) data have compatible dimensions
        if (n != X.getRowCount()) {
            throw new RuntimeException("Illegal vector length.");
        }

        // Solve for the regression coefficients using ordinary least-squares
        double[][] xx = new double[X.getRowCount()][X.getColumnCount() + 1];

        for (int i = 0; i < X.getRowCount(); i++) {

            xx[i][0] = 1.0;
            for (int j = 0; j < X.getColumnCount(); j++) {
                xx[i][j + 1] = X.getData(i, j);
            }
        }

        Matrix XX = new Matrix(xx);
        Matrix XXt = XX.transpose();
        Matrix XTXI = XXt.times(XX);

        Matrix XTXIinv = XTXI.invert();

        // get reciprocal condition of matrix 
        double rcond = 1 / (XTXIinv.norm1() * XTXI.norm1());

        if (inverse == true) {
            if (rcond < this.rcondUpper || Double.isNaN(rcond)) {
                double err = 1e-26;
                return err;
            }
        } else {
            if (rcond < this.rcondUpper || Double.isNaN(rcond)) {
                double err = 100000000;
                return err;
            }
        }

        double[][] p = new double[lengthP][1];
        for (int i = 0; i < lengthP; i++) {
            p[i][0] = modelProp.get(i);
        }

        Matrix PP = new Matrix(p);
        Matrix a = XTXIinv.times(XXt);
        Matrix Coefficients = a.times(PP);
        Matrix PP_hat = XX.times(Coefficients);
        Matrix residuals = PP.minus(PP_hat);

        Matrix b = residuals.transpose().times(residuals);
        double err = Math.sqrt(b.getData(0, 0) / (n - k - 1));

        if (Double.isNaN(err) || Double.isInfinite(err)) {
            logger.info("err value is infinite");
            if (inverse == true) {
                err = 1e-26;
                return err;
            } else {
                err = 100000000;
                return err;
            }
        }
        return err;
    }

    /* rms_inv evaluates selection of descriptors using multiple linear 
     regression analysis. */
    private Double rms_inv(int[] vec) {
        return rms(vec, true);
    }

    /* rmsr replaces all descriptors from Mat (except the ones in the given 
     descriptor vector) in the given position pth and returns the vector with 
     lower S_err, using multiple linear regression analysis. */
    private double[] rmsr(int[] vec, Integer pth) {
        return rmsr(vec, pth, false);
    }

    /* rmsr_inv replaces all descriptors from Mat (except the given descriptor
     vector) in the given position pth and returns the vector with higher s_err, 
     using multiple linear regression analysis. */
    private double[] rmsr_inv(int[] vec, Integer pth) {
        return rmsr(vec, pth, true);
    }

    private double[] rmsr(int[] vec, Integer pth, Boolean inv) {

        Boolean inverse = false;
        if (inv != null) {
            inverse = inv;
        }
        Matrix modelM = this.getModelDescriptor();

        int n_m = modelM.getColumnCount();

        int[] num = new int[n_m];
        int inum = 0;
        List<Integer> Vec = new ArrayList<>();
        for (int i = 0; i < vec.length; i++) {
            Vec.add(vec[i]);
        }

        for (int i = 0; i < n_m; i++) {
            if (Vec.contains(i) && !pth.equals(i)) {
            } else {
                if (this.descrExclude[i] == false) {
                    num[inum] = i;
                    inum = inum + 1;
                }
            }
        }

        int n_n = num.length;
        int desc = 1;
        // A very big number compared to a normal S necessary just to start the program
        double smin = 1000000000;
        double smax = 0.00000000000000000000000000000001;

        for (int j = 0; j < n_n; j++) {
            vec[pth] = num[j];
            double serr = rms(vec);

            if (inverse == true) {
                if (serr > smax) {
                    smax = serr;
                    desc = num[j];
                }
            } else {
                if (serr < smin) {
                    smin = serr;
                    desc = num[j];
                }
            }
        }

        vec[pth] = desc;

        double[] a = new double[vec.length + 1];
        if (inverse == true) {
            a[0] = smax;
        } else {
            a[0] = smin;
        }

        for (int j = 0; j < vec.length; j++) {
            a[j + 1] = vec[j];
        }

        return a;
    }

    public void testMatrix() {

        double[][] m = {{2.0, 0.0, 1.0},
        {-1.0, 1.0, 0.0},
        {-3.0, 3.0, 0.0}};

        Matrix M = new Matrix(m);

        M.show();

        logger.debug("norm-1: \t\t" + M.norm1());
        logger.debug("norm-2: \t\t" + M.norm2());
        logger.debug("frobenius norm: \t" + M.normF());
        logger.debug("normInf: \t\t" + M.normInf());

        Matrix Mt = M.transpose();
        Matrix MM = Mt.times(M);
        Matrix Minv = MM.invert();

        logger.debug("M'*M norm1: " + Minv.norm1());
        logger.debug("inverse M'*M norm1: " + Minv.norm1());

        double err = 1 / (Minv.norm1() * MM.norm1());

        logger.debug("s_err: " + err);
    }

    public void showParameter() {

        if (vecSingle == null) {
            return;
        }
        int[] vec = this.vecSingle;

        System.out.println(ColorStatic.GREEN + "Pexp values" + ColorStatic.RESET);
        System.out.println(this.getModelProperty().toString());

        System.out.println(ColorStatic.GREEN + "Pcalc values" + ColorStatic.RESET);
        System.out.println(this.getModelResponse(this.getModelDescriptor()).toString());

        logger.info("\n --------------------------------- \n"
                + "descriptors selected: " + Arrays.toString(vec) + "\n"
                + "predicted parameters:             \n"
                + "r square                     " + this.modelRSquare() + " \n"
                + "predictive r square          " + this.predictiveRSquare() + " \n");
    }

    public Double RSS() {
        return RSS(this.getModelProperty(), this.getModelResponse(this.getModelDescriptor()));
    }

    public Double RSS(List<Double> pP, List<Double> pp_hat) {

        if (pP.isEmpty()) {
            logger.error("Empty property vector");
            return null;
        }
        double n = pP.size();
        double sumP = 0.0;

        for (Double pP1 : pP) {
            sumP += pP1;
        }
        double meanP = sumP / n;

        double normP = 0.0;
        for (Double P1 : pp_hat) {
            normP += (P1 - meanP) * (P1 - meanP);
        }
        return normP;
    }

    private double rsquare(List<Double> pP, List<Double> pp_hat) {
        double rsq = 0.0;

        double n = pP.size();
        double ySum = 0.0;

        for (Double pP1 : pP) {
            ySum += pP1;
        }

        double yBar = ySum / n;
        double yyBar = 0.0;
        for (Double pP1 : pP) {
            yyBar += (pP1 - yBar) * (pP1 - yBar);
        }

        double rss = 0.0;
        double ssr = 0.0;
        for (int i = 0; i < pP.size(); i++) {
            rss += (pp_hat.get(i) - pP.get(i)) * (pp_hat.get(i) - pP.get(i));
            ssr += (pp_hat.get(i) - yBar) * (pp_hat.get(i) - yBar);
        }
        rsq = ssr / yyBar;
        return rsq;
    }

    public Double TSS() {
        return TSS(this.getModelProperty());
    }

    public Double TSS(List<Double> pP) {

        if (pP.isEmpty()) {
            logger.debug("Empty property vector");
            return null;
        }

        double meanP = 0.0;
        for (Double P1 : pP) {
            meanP += P1;
        }
        meanP = meanP / pP.size();
        double normP = 0.0;
        for (Double P1 : pP) {
            normP = normP + Math.pow(Math.abs(P1 - meanP), 2);
        }
        return normP;
    }

    public Double predictiveRSquare() {

        Double predrsquare = 0.00;

        try {
            List<Double> yPred = this.getModelResponse(this.getPredictDescriptor());
            List<Double> y = this.getPredictProperty();
            Double sumA = 0.0;
            Double sumB = 0.0;
            Double meanY = MathUtil.mean(this.getModelProperty());

            for (int i = 0; i < yPred.size(); i++) {
                sumA = sumA + Math.pow((yPred.get(i) - y.get(i)), 2);
                sumB = sumB + Math.pow((y.get(i) - meanY), 2);
            }
            predrsquare = 1 - (sumA / sumB);
        } catch (Exception e) {
            logger.error(e.getMessage());
        }

        return predrsquare;
    }

    public Double modelRSquare() {

        Double modRSS = null;
        Double modTSS = null;
        Double res = null;

        try {
            modRSS = this.RSS(this.getModelProperty(), this.getModelResponse(this.getModelDescriptor()));
            modTSS = this.TSS(this.getModelProperty());
            res = modRSS / modTSS;
        } catch (Exception e) {
            logger.error(e.getMessage());
        }

        if (modRSS == null || modTSS == null) {
            logger.error("rsquare RSS/TSS could not be calculated");
            return null;
        }
        return res;
    }
}
