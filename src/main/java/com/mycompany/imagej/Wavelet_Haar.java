/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package com.mycompany.imagej;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.io.SaveDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import lib.ImageAccess;

/**
 * A template for processing each pixel of either
 * GRAY8, GRAY16, GRAY32 or COLOR_RGB images.
 *
 * @author Johannes Schindelin
 */
public class Wavelet_Haar implements PlugInFilter {
    int k;                      // Number of nearest neighbors
    int levels;                  // Wavelet decoposition levels
    String referenceImagePath;
    DistanceCalculator distanceCalculator;

    private Map<String, DistanceCalculator> distanceCalculators = new HashMap<>(3);

    public Wavelet_Haar() {
        distanceCalculators.put("L1 (Manhattan)", new L1Calculator());
        distanceCalculators.put("L2 (Euclidean)", new L2Calculator());
        distanceCalculators.put("L-inf (Chebyshev)", new LInfinityCalculator());
    }

    public int setup(String arg, ImagePlus imp) {
        referenceImagePath = imp.getOriginalFileInfo().getFilePath();
        ImageConverter ic = new ImageConverter(imp);
        ic.convertToGray8();
        return DOES_ALL;
    }

    public void run(ImageProcessor img) {

        GenericDialog gd = new GenericDialog("k-nearest neighbor search", IJ.getInstance());
        gd.addNumericField("Number of nearest neighbors (K):", 1, 0);
        gd.addNumericField("Wavelet decomposition level:", 1, 0);
        gd.addChoice("Distance function", new String[] { "L1 (Manhattan)", "L2 (Euclidean)", "L-inf (Chebyshev)"}, "L2 (Euclidean)");
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        k = (int) gd.getNextNumber();
        levels = (int) gd.getNextNumber();
        distanceCalculator = distanceCalculators.get(gd.getNextChoice());

        SaveDialog sd = new SaveDialog("Open search folder...", "any file (required)", "");
        if (sd.getFileName()==null) return;
        String dir = sd.getDirectory();

        search(dir);
    }

    public void search(String dir) {
        IJ.log("");
        IJ.log("Searching images");
        if (!dir.endsWith(File.separator))
            dir += File.separator;
        String[] list = new File(dir).list();  /* lista de arquivos */
        if (list==null) return;

		// Matrix of feature vectors. Each feature vector is an n x 2 (n = number of bands, 2 = number of descriptors) array
        double[][][] featVectors = new double[list.length][3 * levels + 1][2];

        // Array of descriptors of searched images
        Descriptor[] descriptors = new Descriptor[list.length];

        for (int i=0; i<list.length; i++) {
            IJ.showStatus(i+"/"+list.length+": "+list[i]);   /* mostra na interface */
            IJ.showProgress((double)i / list.length);  /* barra de progresso */
            File f = new File(dir+list[i]);

            if (!f.isDirectory()) {
                ImagePlus image = new Opener().openImage(dir, list[i]); /* abre imagem image */

                if (image != null) {
                    // image.show();
                    ImageAccess input = new ImageAccess(image.getProcessor());
                    ImageAccess output = normalize(waveletHaar(input, levels)); // scale to [0, 255]
                    // output.show("Result (Haar Wavelet Transform, " + levels + " levels)");

                    descriptors[i] = new Descriptor(list[i], calculateFeatureVector(output, levels));
                    appendFeatVector(output, levels, featVectors[i]);
                }
            }
        }

        writeFeatVectors(featVectors);

        IJ.showProgress(1.0);

        IJ.log("Calculating distances...");

        // Get reference image
        ImagePlus image = new Opener().openImage(referenceImagePath);
        if (image == null) {
            IJ.log("Image reference not found at path " + referenceImagePath);
            return;
        }

        // Get reference image's descriptor
        ImageAccess input = new ImageAccess(image.getProcessor());
        ImageAccess output = waveletHaar(input, levels);
        Descriptor referenceImageDescriptor = new Descriptor(referenceImagePath, calculateFeatureVector(normalize(output), levels));

        try {
            Files.createDirectories(Paths.get("results"));
        } catch (IOException e) {
            IJ.log(e.getMessage());
            e.printStackTrace();
        }
        IJ.log("Saving results to " + Paths.get("results").toAbsolutePath());

        writeDescriptor(0, referenceImageDescriptor, false, distanceCalculator.getName());

        // Calculate distances to reference descriptor and sort by distance
        SortedMap<Double, Descriptor> map = new TreeMap<>();
        for (Descriptor descriptor : descriptors) {
            double distance = distanceCalculator.calculateDistance(referenceImageDescriptor, descriptor);
            map.put(distance, descriptor);
        }

        // get only k-nearest neighbors
        int i = 0;
        for (Map.Entry<Double, Descriptor> entry : map.entrySet()) {
            writeDescriptor(entry.getKey(), entry.getValue(), true, distanceCalculator.getName());
            showImage(dir, entry.getValue());

            if (++i == k) break;
        }
    }

    static public ImageAccess waveletHaar(ImageAccess input, int levels) {
        int nx = input.getWidth();
        int ny = input.getHeight();
        double[] lowHigh;
        ImageAccess aux1 = new ImageAccess(nx, ny);
        ImageAccess aux2 = new ImageAccess(nx, ny);
        ImageAccess out = input.duplicate();
        ImageAccess temp;

        for (int i = 0; i < levels; i++) {
            int div = 1 << i;

            // vertical split
            for (int y = 0; y < (ny / div); y++) {
                for (int x = 0; x < (nx / div); x += 2) {
                    lowHigh = getLowHigh(out, x, y, 1, 0);
                    aux1.putPixel(x / 2, y, lowHigh[0]);
                    aux1.putPixel(x / 2 + nx / (div << 1), y, lowHigh[1]);
                }
            }

            // horizontal split
            for (int x = 0; x < (nx / div); x++) {
                for (int y = 0; y < (ny / div); y += 2) {
                    lowHigh = getLowHigh(aux1, x, y, 0, 1);
                    aux2.putPixel(x, y / 2, lowHigh[0]);
                    aux2.putPixel(x, y / 2 + ny / (div << 1), lowHigh[1]);
                }
            }

            // copy higher bands
            for (int x = 0; x < nx; x++) {
                for (int y = 0; y < ny; y++) {
					if (x < (nx / div) && y < (ny / div)) continue;

                    aux2.putPixel(x, y, out.getPixel(x, y));
                }
            }

            // rotate references
            temp = out;
            out = aux2;
            aux2 = aux1;
            aux1 = temp;
        }

        return out;
    }

    static final private double SQRT2 = Math.sqrt(2);

    static private double[] getLowHigh(ImageAccess img, int x, int y, int dx, int dy) {
        double first = img.getPixel(x, y);
        double second = img.getPixel(x + dx, y + dy);

        double low = (first + second) / 2 * SQRT2;
        double high = (first - second) / 2 * SQRT2;

        return new double[] { low, high };
    }

    static private ImageAccess normalize(ImageAccess image) {
        return normalize(image, 0, 255);
    }

    static private ImageAccess normalize(ImageAccess image, double newMin, double newMax) {
        double min = image.getMinimum();
        double max = image.getMaximum();
        double div = max - min;
        double fac = newMax - newMin;
        double[][] pixels = image.getArrayPixels();

        for (int i = 0; i < pixels.length; i++) {
            for (int j = 0; j < pixels[i].length; j++) {
                pixels[i][j] = (pixels[i][j] - min) / div * fac + newMin;
            }
        }

        return new ImageAccess(pixels);
    }

    static private double[][] appendFeatVector(ImageAccess input, int levels, double[][] featVector) {
        int nx = input.getWidth();
        int ny = input.getHeight();

        for (int y = 0; y < ny; y++)
            for (int x = 0; x < nx; x++) {
                int band = getBandNumber(x, y, nx, ny, levels);
                double P = input.getPixel(x, y);

                // energy
                featVector[band][0] += Math.pow(P, 2);

                // entropy
                featVector[band][1] -= P > 0 ? P * Math.log(P) : 0;
            }

        return featVector;
    }

    /**
     * Calculate feature vector of a transformed image.
     * 
     * @param input transformed image
     * @param levels how many levels of DWT were made
     * @return array with features of each band, in pairs, i.e. array[2*i] is energy and array[2*i+1] is entropy of
     * i-th band of processed image
     */
    static private double[] calculateFeatureVector(ImageAccess input, int levels) {
        int nx = input.getWidth();
        int ny = input.getHeight();
        double[] featVector = new double[6 * levels + 2];

        for (int y = 0; y < ny; y++) {
            for (int x = 0; x < nx; x++) {
                int band = getBandNumber(x, y, nx, ny, levels);
                double P = input.getPixel(x, y);

                // energy
                featVector[2 * band] += Math.pow(P, 2);

                // entropy
                featVector[2 * band + 1] -= P > 0 ? P * Math.log(P) : 0;
            }
        }

        return featVector;
    }

    static private void showImage(String dir, Descriptor descriptor) {
        new Opener().openImage(dir, descriptor.getId()).show();
    }

    /**
	 * Gets the band number based on the coordinates of a pixel and the number of decomposition level, labeled in read-order, for each level of decomposition.
     * Band labels example for an image with 3 levels of decomposition:
     * -------------------------
	 * |0 |1 |  4  |           |
	 * |2 |3 |     |           |
	 * ------------|     7     |
	 * |  5  |  6  |           |
	 * |     |     |           |
     * -------------------------
	 * |		   |		   |
	 * |		   |		   |
	 * |	 8	   |	 9	   |
	 * |		   |		   |
	 * |		   |		   |
     * -------------------------
     */
    static private int getBandNumber(int x, int y, int nx, int ny, int levels) {
        int localBand = 0;
        for (int level = 0; level < levels; level++) {
            int divX = nx / 2;
            int divY = ny / 2;

            if (x < divX && y < divY) {
                localBand = 0; // LL band
            } else if (x >= divX && y < divY) {
                localBand = 1; // LH band
                x -= divX;
            } else if (x < divX && y >= divY) {
                localBand = 2; // HL band
                y -= divY;
            } else {
                localBand = 3; // HH band
                x -= divX;
                y -= divY;
            }

            // half the size
            nx = divX;
            ny = divY;

            // if the pixel is not in the LL band, we do not decompose further
            if (localBand != 0) {
                return (levels - (level + 1)) * 3 + localBand;
            }
        }

        return localBand;
    }

    /**
     * Write descriptor's info to a file and to ImageJ's GUI.
     * 
     * @param distance distance of this image to reference image
     * @param descriptor Descriptor object of this image
     * @param append whether to append to end of file or not
     */
    static private void writeDescriptor(double distance, Descriptor descriptor, boolean append, String distanceFunction) {
        try {
            FileWriter writer = new FileWriter("results/knn.txt", append);
            StringBuilder s = new StringBuilder();
            double[] featureVector = descriptor.getFeatureVector();

            if (!append) {
                s.append("Distance function: " + distanceFunction + "\n");
            }

            s.append("Image: " + descriptor.getId() + ", distance: " + distance + ", features vector: (");

            for (int i = 0; i < featureVector.length / 2; i++) {
                double energy = featureVector[2 * i];
                double entropy = featureVector[2 * i + 1];
                s.append(String.format("%.0f", energy) + ", " + String.format("%.0f", entropy) + ", ");
            }
            s.deleteCharAt(s.length() - 1);
            s.deleteCharAt(s.length() - 1);
            s.append(")\n");

            String str = s.toString();
            writer.write(str);
            IJ.log(str);

            writer.close();
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    static private void writeFeatVectors(double[][][] featVectors) {
        try {
            FileWriter myWriter = new FileWriter("results/feat_vectors.txt");
            String s = "";
            for (double[][] featVector : featVectors) {
                for (double[] energy_entropy : featVector) {
                    double energy = energy_entropy[0];
                    double entropy = energy_entropy[1];
                    s += "(" + String.format("%.0f", energy) + "," + String.format("%.0f", entropy) + ")" + ", ";
                }
                s += "\n";
            }
            myWriter.write(s);
            myWriter.close();
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    /**
     * Main method for debugging.
     *
     * For debugging, it is convenient to have a method that starts ImageJ, loads
     * an image and calls the plugin, e.g. after setting breakpoints.
     *
     * @param args unused
     */
    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = Wavelet_Haar.class;
        java.net.URL url = clazz.getProtectionDomain().getCodeSource().getLocation();
        java.io.File file = new java.io.File(url.toURI());
        System.setProperty("plugins.dir", file.getAbsolutePath());

        // start ImageJ
        new ImageJ();

        // open the sample
		OpenDialog sd = new OpenDialog("Open image...");
		if (sd.getFileName()==null) return;
        ImagePlus image = IJ.openImage(sd.getPath());
        image.show();

        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }

    /**
     * Class that describes a certain image with an id and featureVector.
     */
    private static class Descriptor {
        private String _id;
        private double[] _featureVector;

        public Descriptor(String id, double[] featureVector) {
            _id = id;
            _featureVector = featureVector;
        }

        public String getId() {
            return _id;
        }

        /**
         * Getter for featureVector
         * @return a copy of internal featureVector
         */
        public double[] getFeatureVector() {
            return _featureVector.clone();
        }
    }

    private interface DistanceCalculator {
        public double calculateDistance(Descriptor a, Descriptor b);
        public String getName();
    }

    private class L1Calculator implements DistanceCalculator {

        @Override
        public double calculateDistance(Descriptor a, Descriptor b) {
            double sum = 0;

            for (int i = 0; i < a._featureVector.length; i++) {
                sum += Math.abs(a._featureVector[i] - b._featureVector[i]);
            }

            return sum;
        }

        @Override
        public String getName() {
            return "L1 (Manhattan)";
        }
    }

    private class L2Calculator implements DistanceCalculator {

        @Override
        public double calculateDistance(Descriptor a, Descriptor b) {
            double sum = 0;

            for (int i = 0; i < a._featureVector.length; i++) {
                sum += Math.pow(a._featureVector[i] - b._featureVector[i], 2);
            }

            return Math.sqrt(sum);
        }

        @Override
        public String getName() {
            return "L2 (Euclidean)";
        }
    }

    private class LInfinityCalculator implements DistanceCalculator {

        @Override
        public double calculateDistance(Descriptor a, Descriptor b) {
            double max = 0;

            for (int i = 0; i < a._featureVector.length; i++) {
                max = Math.max(max, Math.abs(a._featureVector[i] - b._featureVector[i]));
            }

            return max;
        }

        @Override
        public String getName() {
            return "L-inf (Chebyshev)";
        }
    }
}
