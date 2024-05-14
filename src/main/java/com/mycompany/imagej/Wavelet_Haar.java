/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package com.mycompany.imagej;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import lib.ImageAccess;

/**
 * A template for processing each pixel of either
 * GRAY8, GRAY16, GRAY32 or COLOR_RGB images.
 *
 * @author Johannes Schindelin
 */
public class Wavelet_Haar implements PlugInFilter {
	protected ImagePlus image;

	// plugin parameters
	public double value;
	public String name;

	@Override
	public int setup(String arg, ImagePlus imp) {
		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}

		image = imp;
		return DOES_8G | DOES_16 | DOES_32 | DOES_RGB;
	}

	@Override
	public void run(ImageProcessor ip) {
		if (showDialog()) {
			ImageAccess out = process(new ImageAccess(ip));
			out.show("Result");;
		}
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Wavelet Haar");

		gd.addNumericField("levels", 1, 0);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		value = gd.getNextNumber();

		return true;
	}

	static public ImageAccess process(ImageAccess input) {
		return waveletHaar(input, 1);
	}

	static public ImageAccess waveletHaar(ImageAccess input, int max) {
		int nx = input.getWidth();
		int ny = input.getHeight();
		double[] lowHigh;
		ImageAccess middle = new ImageAccess(nx, ny);
		ImageAccess out = new ImageAccess(nx, ny);

		// Low Pass | High Pass
		for (int x = 0; x < nx - 1; x += 2) {
			for (int y = 0; y < ny; y++) {
				lowHigh = getLowHigh(input, x, y, 1, 0);
				middle.putPixel(x / 2, y, lowHigh[0]);
				middle.putPixel((x + nx) / 2, y, lowHigh[1]);
			}
		}

		// LL HL
		// LH HH
		for (int x = 0; x < nx; x++) {
			for (int y = 0; y < ny - 1; y +=2) {
				lowHigh = getLowHigh(middle, x, y, 0, 1);
				out.putPixel(x, y / 2, lowHigh[0]);
				out.putPixel(x, (y + ny) / 2, lowHigh[1]);
			}
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

	public void showAbout() {
		IJ.showMessage("WaveletHaar",
				"a plugin to compute the Wavelet Haar descriptors of an image");
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

		// open the Clown sample
		ImagePlus image = IJ.openImage("http://imagej.net/images/clown.jpg");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
