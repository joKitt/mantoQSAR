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
import java.util.List;
import java.util.Random;

public abstract class Vector {

    public static Double norm(Double[] vector) {

        Double a = 0.0;
        for (Double vector1 : vector) {
            a = a + (Math.pow(vector1, 2));
        }

        return Math.sqrt(a);
    }

    public static Double distance(Double[] vec1, Double[] vec2) {

        Double n1 = Vector.norm(vec1);
        Double n2 = Vector.norm(vec2);

        for (int i = 0; i < 3; i++) {
            vec1[i] = vec1[i] / n1;
            vec2[i] = vec2[i] / n2;
        }

        Double a = 0.0;
        for (int i = 0; i < 3; i++) {
            a = a + Math.pow((vec1[i] - vec2[i]), 2);
        }
        return Math.sqrt(a);
    }

    public static Double[] crossProduct(Double[] vec1, Double[] vec2) {

        Double[] vecOut = new Double[3];

        vecOut[0] = vec1[1] * vec2[2] - vec2[1] * vec1[2];
        vecOut[1] = vec1[2] * vec2[0] - vec2[2] * vec1[0];
        vecOut[2] = vec1[0] * vec2[1] - vec2[0] * vec1[1];

        return vecOut;
    }

    public static double cosineSimilarity(double[] vectorA, double[] vectorB) {
        double dotProduct = 0.0;
        double normA = 0.0;
        double normB = 0.0;
        for (int i = 0; i < vectorA.length; i++) {
            dotProduct += vectorA[i] * vectorB[i];
            normA += Math.pow(vectorA[i], 2);
            normB += Math.pow(vectorB[i], 2);
        }
        return dotProduct / (Math.sqrt(normA) * Math.sqrt(normB));
    }

    public static List<Double[]> getSphere120() {

        List<Double[]> sphere = new ArrayList<>();

        sphere.add(new Double[]{0.173274508128383, -0.931181673003991, 0.320743880213411});
        sphere.add(new Double[]{0.852248111355265, -0.0160141095450480, -0.522892632370034});
        sphere.add(new Double[]{0.450672918920059, -0.870494468471259, 0.197821385378360});
        sphere.add(new Double[]{-0.763577001293199, 0.615366152063854, 0.195639111606579});
        sphere.add(new Double[]{0.994288203446936, 0.105025035187744, 0.0189923792632013});
        sphere.add(new Double[]{0.0297983146546386, 0.177157390519868, 0.983731324817876});
        sphere.add(new Double[]{-0.392598542801034, 0.865523644414814, -0.311022837022307});
        sphere.add(new Double[]{-0.103106604621516, -0.0338676688588762, -0.994093561537086});
        sphere.add(new Double[]{0.936907863799216, 0.234611495189493, -0.259154589147367});
        sphere.add(new Double[]{0.159925349128495, 0.906067771125104, 0.391746189304395});
        sphere.add(new Double[]{0.244890464307378, 0.941661438063290, -0.230873117871904});
        sphere.add(new Double[]{0.393718190898195, 0.0146911488217922, -0.919113788549681});
        sphere.add(new Double[]{-0.561974396781707, 0.826423439144713, -0.0347717787014389});
        sphere.add(new Double[]{0.372683587582338, -0.829301868990737, -0.416371653256113});
        sphere.add(new Double[]{0.648849868890085, -0.0779161103683320, -0.756916724208410});
        sphere.add(new Double[]{-0.158306208337498, -0.974455426752739, -0.159297726518433});
        sphere.add(new Double[]{-0.626407959504375, -0.240349284588997, -0.741515535688311});
        sphere.add(new Double[]{0.971801945124045, -0.194488514879477, 0.133323655189534});
        sphere.add(new Double[]{-0.148697328071566, 0.914869727558839, 0.375369799292368});
        sphere.add(new Double[]{-0.901918287999188, -0.181250602441069, -0.392035229140651});
        sphere.add(new Double[]{-0.217683432827486, -0.659363767544227, -0.719620278426285});
        sphere.add(new Double[]{0.354342991954553, 0.391218035759914, -0.849346509116774});
        sphere.add(new Double[]{0.601966680138635, 0.259477796153986, -0.755186989629682});
        sphere.add(new Double[]{-0.0389555391253100, 0.564779275907358, -0.824322045972942});
        sphere.add(new Double[]{0.457332897344287, -0.880809156039010, -0.122563663638647});
        sphere.add(new Double[]{0.845569257033841, -0.447459362327256, 0.291191948077078});
        sphere.add(new Double[]{-0.385824854278436, -0.0891644936209750, -0.918253164926938});
        sphere.add(new Double[]{-0.0558619497280788, 0.692554734177063, 0.719199125932113});
        sphere.add(new Double[]{0.711520631402334, -0.682050220473267, 0.168955283555123});
        sphere.add(new Double[]{0.0349940068154628, 0.827680158439822, -0.560108002810200});
        sphere.add(new Double[]{0.915979777032312, 0.363190709583407, 0.170509696322932});
        sphere.add(new Double[]{-0.958808007964875, -0.276385367047786, -0.0655616713048766});
        sphere.add(new Double[]{0.511368387165166, 0.247282238708304, 0.823015107411466});
        sphere.add(new Double[]{0.340164674284071, -0.383158333303221, -0.858765210048421});
        sphere.add(new Double[]{-0.508993373687831, 0.182376575301183, -0.841227989503020});
        sphere.add(new Double[]{-0.863827966988723, 0.352035773873510, 0.360377659353918});
        sphere.add(new Double[]{0.267941446492095, 0.632914085637931, -0.726379474828980});
        sphere.add(new Double[]{0.577857426602841, 0.814996752390490, 0.0431403304680882});
        sphere.add(new Double[]{0.0753411110006286, -0.743280075197619, -0.664724338961206});
        sphere.add(new Double[]{-0.918991595241381, -0.251284691502909, 0.303826351213963});
        sphere.add(new Double[]{0.858771215416117, -0.176333288567449, 0.481059841304296});
        sphere.add(new Double[]{-0.764127873054394, 0.629727639365092, -0.139827371573081});
        sphere.add(new Double[]{0.0858597020006137, -0.927490239209515, -0.363854322144787});
        sphere.add(new Double[]{0.220840800923263, 0.445339482189347, 0.867699306356122});
        sphere.add(new Double[]{-0.00633965289717515, -0.452700996200572, -0.891639847045966});
        sphere.add(new Double[]{0.589841766452383, 0.765711795498867, -0.256460789950889});
        sphere.add(new Double[]{0.722499131378577, 0.530718586795090, 0.443094557388629});
        sphere.add(new Double[]{-0.858193947954927, 0.0521032285754242, 0.510674457228432});
        sphere.add(new Double[]{0.617818075181250, -0.375886521726449, -0.690659213189636});
        sphere.add(new Double[]{0.828188179949387, -0.353818063352775, -0.434645967009259});
        sphere.add(new Double[]{-0.761291334715235, -0.273147904531209, 0.588069490738728});
        sphere.add(new Double[]{-0.852000160725711, -0.503218451229108, 0.144453855836181});
        sphere.add(new Double[]{0.676181650998956, -0.422481188997310, 0.603562772042598});
        sphere.add(new Double[]{0.115192953915451, 0.234127070225533, -0.965357497694943});
        sphere.add(new Double[]{0.320237629958902, 0.941749816245665, 0.102738230272678});
        sphere.add(new Double[]{-0.330903548145641, 0.745637962194213, 0.578382979659228});
        sphere.add(new Double[]{-0.394298135656359, 0.515354272423213, -0.760880380948977});
        sphere.add(new Double[]{0.622427014584037, 0.574622569482561, -0.531407107740631});
        sphere.add(new Double[]{-0.976542510471729, 0.0598198915966279, 0.206848509327334});
        sphere.add(new Double[]{0.922273979707909, 0.109479228209016, 0.370708787250173});
        sphere.add(new Double[]{-0.525191198781830, -0.780576970864677, -0.338930372315361});
        sphere.add(new Double[]{-0.239768333279278, -0.855963833846581, -0.458079754522224});
        sphere.add(new Double[]{0.489659832795081, 0.546714234297556, 0.679217781100827});
        sphere.add(new Double[]{0.775511706084643, 0.618176740605621, 0.128215097004595});
        sphere.add(new Double[]{0.693187907693943, -0.0399055998947334, 0.719651351505668});
        sphere.add(new Double[]{0.305671776042289, 0.0478362709306130, 0.950934517469319});
        sphere.add(new Double[]{0.546617315668246, -0.691267884740045, 0.472607894283057});
        sphere.add(new Double[]{-0.670364003997458, 0.301624639286742, 0.677963626693674});
        sphere.add(new Double[]{0.494305089446221, 0.775419042627637, 0.392922112991878});
        sphere.add(new Double[]{-0.0679140503795539, -0.844893542597665, 0.530605864496248});
        sphere.add(new Double[]{-0.665413456791141, -0.744683565833608, -0.0516847975579591});
        sphere.add(new Double[]{0.256979238485514, -0.759827020392880, 0.597180517153970});
        sphere.add(new Double[]{-0.119011643510571, -0.0849276776567276, 0.989254021107192});
        sphere.add(new Double[]{-0.810697424932962, -0.537008227647727, -0.233220600817578});
        sphere.add(new Double[]{-0.586293894484417, 0.660858270187243, -0.468535821485850});
        sphere.add(new Double[]{0.478412628343703, -0.209009379712853, 0.852898842907715});
        sphere.add(new Double[]{-0.486292440212382, 0.823042872409809, 0.293462251009729});
        sphere.add(new Double[]{0.0612052574341384, -0.598047043128879, 0.799120548269925});
        sphere.add(new Double[]{-0.419707627837192, -0.166383775912758, 0.892279074196123});
        sphere.add(new Double[]{-0.514961167289599, -0.599808774193097, -0.612408712041809});
        sphere.add(new Double[]{-0.994245080828553, 0.0409714847805815, -0.0989851336469236});
        sphere.add(new Double[]{-0.933641135179573, 0.358195459658722, 0.00320052815553610});
        sphere.add(new Double[]{-0.138571597581445, -0.970747355771586, 0.196080298873190});
        sphere.add(new Double[]{0.156023671555975, -0.987554528631247, -0.0198158243358561});
        sphere.add(new Double[]{0.965026228357058, -0.127142451972627, -0.229257879884944});
        sphere.add(new Double[]{-0.486913383332077, -0.444108725682401, 0.752118871525548});
        sphere.add(new Double[]{-0.639649273954244, -0.721185899579965, 0.265969367741560});
        sphere.add(new Double[]{-0.126965925057487, 0.428701156259558, 0.894480280663590});
        sphere.add(new Double[]{0.368925553822161, 0.792335543377632, -0.485899498288983});
        sphere.add(new Double[]{0.363615042265533, -0.634754233015295, -0.681814611685170});
        sphere.add(new Double[]{0.0181419197341694, 0.998320420157874, 0.0550200821897333});
        sphere.add(new Double[]{-0.416640186659971, 0.486760258279699, 0.767773017121222});
        sphere.add(new Double[]{0.158622176591530, -0.277802496871058, 0.947451728493632});
        sphere.add(new Double[]{-0.165840441811879, -0.412093533084181, 0.895921797842779});
        sphere.add(new Double[]{-0.829774691562244, 0.432617334103461, -0.352585030135990});
        sphere.add(new Double[]{-0.383201978768780, -0.837774375834036, 0.388960587545297});
        sphere.add(new Double[]{0.828140797510811, 0.525820865432354, -0.194152612586409});
        sphere.add(new Double[]{0.637983627109885, -0.624528069902724, -0.450490378857634});
        sphere.add(new Double[]{-0.0838466102201188, 0.955705269626340, -0.282129728251103});
        sphere.add(new Double[]{-0.758228311082698, 0.0297640162063160, -0.651309397761117});
        sphere.add(new Double[]{-0.676530174401932, 0.370581868273105, -0.636377248203380});
        sphere.add(new Double[]{-0.332168898905184, -0.388694773841705, -0.859406885816179});
        sphere.add(new Double[]{-0.289527996094421, -0.677803293486670, 0.675837432239568});
        sphere.add(new Double[]{0.791545469942480, 0.319429817000991, -0.520980192545025});
        sphere.add(new Double[]{-0.281440094376153, 0.958895673581540, 0.0362016638312103});
        sphere.add(new Double[]{-0.646759293418034, -0.00883164607668725, 0.762643047831016});
        sphere.add(new Double[]{-0.654868658544307, -0.557621535109242, 0.510103189206648});
        sphere.add(new Double[]{0.149587222260931, -0.175144250376418, -0.973112611415730});
        sphere.add(new Double[]{-0.371859665426661, 0.177301802776219, 0.911199462226067});
        sphere.add(new Double[]{-0.265819111885755, 0.753762023812495, -0.600985034101815});
        sphere.add(new Double[]{-0.911520146490951, 0.147347303190969, -0.383952854896878});
        sphere.add(new Double[]{0.394039355417787, -0.503654107949121, 0.768807860214702});
        sphere.add(new Double[]{-0.746032032734772, -0.425792038313379, -0.511993502148770});
        sphere.add(new Double[]{-0.203230840851441, 0.296157389813980, -0.933267392436587});
        sphere.add(new Double[]{-0.636195001833723, 0.584122541174077, 0.504040451287521});
        sphere.add(new Double[]{0.897023256848971, -0.431716404937229, -0.0947112579376906});
        sphere.add(new Double[]{0.231848438866377, 0.731984209346307, 0.640660142870529});
        sphere.add(new Double[]{-0.409698680003785, -0.912190143854895, 0.00749219975324191});
        sphere.add(new Double[]{0.756766114363684, 0.262761645772634, 0.598549384480302});
        sphere.add(new Double[]{0.717092059928110, -0.680017690918165, -0.152823158017321});

        return sphere;
    }

    public static List<Double[]> calcSphere(int nPoint) {

        if (nPoint == 120) {
            return Vector.getSphere120();
        }

        Random random = new Random();
        Double k = 0.1; // stepsize
        Double smudge = 0.001;
        Double tol = 0.001; // tolerance setting

        List<Double[]> sphere = new ArrayList<>();
        System.out.println("Sphere of " + nPoint + " points");
        for (int i = 0; i < nPoint; i++) {

            Double[] vec = new Double[3];
            for (int ii = 0; ii < 3; ii++) {

                vec[ii] = random.nextGaussian();
                System.out.print(vec[ii] + "; ");
            }
            sphere.add(vec);
            System.out.println("");
        }

        Double delta = 100.0;

        while (delta > tol) {
            delta = 0.0;

            for (int kk = 0; kk < sphere.size(); kk++) {

                Double[] vec = Arrays.copyOf(sphere.get(kk), 3);

                List<Double[]> dx = new ArrayList<>();
                List<Double> d = new ArrayList<>();

                for (Double[] sphere1 : sphere) {
                    Double[] b = new Double[3];
                    Double a = 0.0;
                    for (int j = 0; j < 3; j++) {
                        b[j] = vec[j] - sphere1[j];
                        a = a + (Math.pow(b[j], 2));
                    }
                    dx.add(b);
                    d.add(a);
                }

                Double fxx = 0.0;
                Double fxy = 0.0;
                Double fxz = 0.0;

                for (int i = 0; i < sphere.size(); i++) {
                    fxx = fxx + (dx.get(i)[0] / (smudge + d.get(i)));
                    fxy = fxy + (dx.get(i)[1] / (smudge + d.get(i)));
                    fxz = fxz + (dx.get(i)[2] / (smudge + d.get(i)));
                }

                Double[] xold = Arrays.copyOf(vec, 3);

                vec[0] = vec[0] + k * fxx;
                vec[1] = vec[1] + k * fxy;
                vec[2] = vec[2] + k * fxz;

                Double nv = Vector.norm(vec);
                vec[0] = vec[0] / nv;
                vec[1] = vec[1] / nv;
                vec[2] = vec[2] / nv;

                Double[] dVec = new Double[3];
                dVec[0] = xold[0] - vec[0];
                dVec[1] = xold[1] - vec[1];
                dVec[2] = xold[2] - vec[2];

                sphere.set(kk, vec);
                delta = delta + (Vector.norm(dVec) / sphere.size());

            }
        } 

        return sphere;
    }

}
