"""
This file is part of the Duh project.
Copyright 2014 David W. Hogg (NYU) and Dustin Lang (CMU).

## bugs
* No output.
* Needs some kind of robust sigma estimation.
* Needs some kind of IRLS or sigma-clipping.
"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
import numpy as np
import astropy.io.fits as pf
import pylab as plt

def populate_A_matrix(ones, image, nhalf=1):
    """
    ## bugs:
    - BROKEN
    - Needs more information in this comment header.
    - Not yet written.
    - There must be a MUCH faster way to write this.
    """
    ny, nx = image.shape
    npix = nx * ny
    nfull = 2 * nhalf + 1
    A = np.zeros((1 + nfull ** 2, npix))
    A[0] = ones.reshape(npix)
    kk = 1
    for ii in range(-nhalf, nhalf + 1):
        y1 = 0 - ii
        y2 = ny - ii
        y3 = 0
        y4 = ny
        if y1 < 0:
            y3 -= y1
            y1 = 0
        if y2 > ny:
            y4 -= (y2 - ny)
            y2 = ny
        for jj in range(-nhalf, nhalf + 1):
            x1 = 0 - jj
            x2 = nx - jj
            x3 = 0
            x4 = nx
            if x1 < 0:
                x3 -= x1
                x1 = 0
            if x2 > nx:
                x4 -= (x2 - nx)
                x2 = nx
            foo = np.zeros((ny, nx))
            foo[y1: y2, x1: x2] = image[y3: y4, x3: x4]
            A[kk] = foo.reshape(npix)
            kk += 1
    return A

def synthesize_image(image_to_synthesize, weight_image, basis_image):
    """
    ## input:
    - `image_to_synthesize` - new image in which we are looking for anomalies
    - `weight_image` - image of inverse variance values or 1 (to use in fit) and 0 (to not use) values.
    - `basis_image` - old image that represents "truth" or whatever

    ## output:
    - `synthetic_image` - best-fit image created from `basis_image` to approximate `image_to_synthesize`

    ## comments:
    - Images must be input with the same size and shape and very close to perfectly aligned
      (this is really a "bug" not a "comment")!

    ## bugs:
    - Not yet written.
    - Ought to have user-settable order and options related to model complexity.
    - Can't take multiple `basis_image`s.  That would be cool!
    """
    assert weight_image.shape == image_to_synthesize.shape
    assert basis_image.shape == image_to_synthesize.shape
    ny, nx = image_to_synthesize.shape
    xg = (np.arange(nx) + 0.5 - 0.5 * nx) * 2. / nx
    yg = (np.arange(ny) + 0.5 - 0.5 * ny) * 2. / ny
    yg, xg = np.meshgrid(xg, yg)
    A = np.vstack((populate_A_matrix(np.ones((ny, nx)), basis_image, nhalf=3),
                   populate_A_matrix(xg, xg * basis_image),
                   populate_A_matrix(yg, yg * basis_image)))
    npix = ny * nx
    W = weight_image.reshape(npix)
    ATA = np.dot(W[None, :] * A, A.T)
    ATb = np.dot(W[None, :] * A, image_to_synthesize.reshape(npix))
    pars = np.linalg.solve(ATA, ATb)
    return np.dot(A.T, pars).reshape((ny, nx))

def iterative_reweight_synthesize_image(data_j2c, invvar, data_j3c):
    """
    ## bugs:
    - DOES NOTHING.
    - Total hack.
    """
    new_invvar = invvar
    for ii in range(5):
        diff = synthesize_image(data_j2c, new_invvar, data_j3c) - data_j2c # use new invvar here
        chi2 = invvar * (diff ** 2) # use original invvar here
        factor = (chi2 + 25.) / 25. # MAGIC 25.
        new_invvar = invvar / factor # update new invvar here
        print "median invvar, chi2, new_invvar", np.median(invvar), np.median(chi2), np.median(new_invvar)
    return synthesize_image(data_j2c, new_invvar, data_j3c), new_invvar

def synthesize_and_plot(fn1, fn2, fn):
    data_j2c = pf.open(fn1)[0].data
    data_j3c = pf.open(fn2)[0].data

    sigma_j2c = np.median(data_j2c) - np.percentile(data_j2c, 16.)
    invvar = np.zeros_like(data_j2c) + 1. / (sigma_j2c * sigma_j2c)
    synth_j2c_j3c, new_invvar = iterative_reweight_synthesize_image(data_j2c, invvar, data_j3c)
    resid_j2c_j3c = data_j2c - synth_j2c_j3c

    plt.clf()
    plt.subplot(231)
    plt.gray()
    vmin_j2c = np.percentile(data_j2c, 1.)
    vmax_j2c = np.percentile(data_j2c, 99.)
    y1 = 200
    y2 = y1 + 400
    x1 = 200
    x2 = x1 + 400
    plt.imshow(data_j2c[y1: y2, x1: x2], interpolation="nearest", vmin=vmin_j2c, vmax=vmax_j2c)
    plt.title(fn1)

    plt.subplot(232)
    vmin_j3c = np.percentile(data_j3c, 5.)
    vmax_j3c = np.percentile(data_j3c, 95.)
    plt.imshow(data_j3c[y1: y2, x1: x2], interpolation="nearest", vmin=vmin_j3c, vmax=vmax_j3c)
    plt.title(fn2)

    plt.subplot(233)
    plt.imshow(np.sqrt(new_invvar)[y1: y2, x1: x2], interpolation="nearest")
    plt.title("sqrt inverse variance")

    plt.subplot(234)
    plt.imshow(synth_j2c_j3c[y1: y2, x1: x2], interpolation="nearest", vmin=vmin_j2c, vmax=vmax_j2c)
    plt.title("synthetic " + fn1)

    plt.subplot(235)
    diff = np.median(data_j2c) - np.median(resid_j2c_j3c)
    plt.imshow(resid_j2c_j3c[y1: y2, x1: x2], interpolation="nearest", vmin=vmin_j2c - diff, vmax=vmax_j2c - diff)
    plt.title("residual")

    plt.subplot(236)
    plt.imshow((np.sqrt(new_invvar) * resid_j2c_j3c)[y1: y2, x1: x2], interpolation="nearest", vmin=-5., vmax=5.)
    plt.title("chi image")

    plt.savefig(fn)
    return None

if __name__ == "__main__":
    synthesize_and_plot("j2c.fits", "j3c.fits", "compare_j2c_j3c.png")
    synthesize_and_plot("j3c.fits", "j2c.fits", "compare_j3c_j2c.png")

if False:
    new = np.random.normal(size=(1700, 1500))
    weight = np.ones((1700, 1500))
    old = np.random.uniform(size=(1700, 1500))
    print new - synthesize_image(new, weight, old)
