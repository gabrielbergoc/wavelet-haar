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

    public int setup(String arg, ImagePlus imp) {
        ImageConverter ic = new ImageConverter(imp);
        ic.convertToGray8();
        return DOES_ALL;
    }

    public void run(ImageProcessor img) {

        GenericDialog gd = new GenericDialog("k-nearest neighbor search", IJ.getInstance());
        gd.addNumericField("Number of nearest neighbors (K):", 1, 0);
        gd.addNumericField("Wavelet decomposition level:", 1, 0);
        gd.showDialog();
        if (gd.wasCanceled())
            return;
        k = (int) gd.getNextNumber();
        levels = (int) gd.getNextNumber();

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

        for (int i=0; i<list.length; i++) {
            IJ.showStatus(i+"/"+list.length+": "+list[i]);   /* mostra na interface */
            IJ.showProgress((double)i / list.length);  /* barra de progresso */
            File f = new File(dir+list[i]);
            if (!f.isDirectory()) {
                ImagePlus image = new Opener().openImage(dir, list[i]); /* abre imagem image */
                if (image != null) {
					// image.show();
                    ImageAccess input = new ImageAccess(image.getProcessor());
					ImageAccess output = waveletHaar(input, levels);
					// output.show("Result (Haar Wavelet Transform, " + levels + " levels)");

					appendFeatVector(output, levels, featVectors[i]);
                }
            }
        }

		writeFeatVectors(featVectors);

        IJ.showProgress(1.0);
        IJ.showStatus("");      
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

	static private double[][] appendFeatVector(ImageAccess input, int levels, double[][] featVector) {
		int nx = input.getWidth();
		int ny = input.getHeight();

		for (int y = 0; y < ny; y++)
		for (int x = 0; x < nx; x++) {
			int band = getBandNumber(x, y, nx, ny, levels);
			double P = input.getPixel(x, y);
			P += 128; // change to normalize with max and min

			// energy
			featVector[band][0] += Math.pow(P, 2);

			// entropy
			featVector[band][1] -= P * Math.log(P);
		}

		return featVector;
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
}
