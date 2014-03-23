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
    A = np.vstack((populate_A_matrix(np.ones((ny, nx)), basis_image, nhalf=6),
                   populate_A_matrix(xg, xg * basis_image, nhalf=1),
                   populate_A_matrix(yg, yg * basis_image, nhalf=1)))
    npix = ny * nx
    W = weight_image.reshape(npix)
    ATA = np.dot(W[None, :] * A, A.T)
    ATb = np.dot(W[None, :] * A, image_to_synthesize.reshape(npix))
    pars = np.linalg.solve(ATA, ATb)
    return np.dot(A.T, pars).reshape((ny, nx))

def iterative_reweight_synthesize_image(data1, invvar, data2):
    """
    ## bugs:
    - DOES NOTHING.
    - Total hack.
    """
    new_invvar = invvar
    for ii in range(5):
        diff = synthesize_image(data1, new_invvar, data2) - data1 # use new invvar here
        chi2 = invvar * (diff ** 2) # use original invvar here
        factor = (chi2 + 25.) / 25. # MAGIC 25.
        new_invvar = invvar / factor # update new invvar here
        I = invvar > 0
        print "median invvar, chi2, new_invvar", np.median(invvar[I]), np.median(chi2[I]), np.median(new_invvar[I])
    return synthesize_image(data1, new_invvar, data2), new_invvar

def synthesize_and_plot(fn1, fn2, fn, maskzero=False):
    """
    ## bugs:
    - No comment header.
    """
    data1 = pf.open(fn1)[0].data
    data2 = pf.open(fn2)[0].data

    sigma1 = np.median(data1[data1 != 0]) - np.percentile(data1[data1 != 0], 16.)
    invvar = np.zeros_like(data1) + 1. / (sigma1 * sigma1)
    invvar[~np.isfinite(data1)] = 0.
    data1[~np.isfinite(data1)] = 0.
    data2[~np.isfinite(data2)] = 0.
    if maskzero:
        invvar[data1 == 0] = 0.
        invvar[data2 == 0] = 0.
    synth12, new_invvar = iterative_reweight_synthesize_image(data1, invvar, data2)
    resid12 = data1 - synth12

    plt.figure(figsize=(9, 6), dpi=600)
    plt.clf()
    plt.subplot(231)
    plt.gray()
    vmin1 = np.percentile(data1, 1.)
    vmax1 = np.percentile(data1, 99.)
    ny, nx = data1.shape
    y1 = 0
    y2 = ny
    x1 = 0
    x2 = nx
    plt.imshow(data1[y1: y2, x1: x2], interpolation="nearest", vmin=vmin1, vmax=vmax1)
    plt.title(fn1)

    plt.subplot(232)
    vmin2 = np.percentile(data2, 5.)
    vmax2 = np.percentile(data2, 95.)
    plt.imshow(data2[y1: y2, x1: x2], interpolation="nearest", vmin=vmin2, vmax=vmax2)
    plt.title(fn2)

    plt.subplot(233)
    plt.imshow(np.sqrt(new_invvar)[y1: y2, x1: x2], interpolation="nearest")
    plt.title("sqrt inverse variance")

    plt.subplot(234)
    plt.imshow(synth12[y1: y2, x1: x2], interpolation="nearest", vmin=vmin1, vmax=vmax1)
    plt.title("synthetic " + fn1)

    plt.subplot(235)
    diff = np.median(data1) - np.median(resid12)
    plt.imshow(resid12[y1: y2, x1: x2], interpolation="nearest", vmin=vmin1 - diff, vmax=vmax1 - diff)
    plt.title("residual")

    plt.subplot(236)
    plt.imshow((np.sqrt(new_invvar) * resid12)[y1: y2, x1: x2], interpolation="nearest", vmin=-5., vmax=5.)
    plt.title("chi image")

    plt.savefig(fn)
    return None

if __name__ == "__main__":
    synthesize_and_plot("j2c.fits", "j3c.fits", "compare_j2c_j3c.pdf")
    synthesize_and_plot("j3c.fits", "j2c.fits", "compare_j3c_j2c.pdf")

if False:
    synthesize_and_plot("h.fits", "j2c.fits", "compare_h_j2c.png", maskzero=True)
    synthesize_and_plot("h.fits", "j3c.fits", "compare_h_j3c.png", maskzero=True)

if False:
    new = np.random.normal(size=(1700, 1500))
    weight = np.ones((1700, 1500))
    old = np.random.uniform(size=(1700, 1500))
    print new - synthesize_image(new, weight, old)
