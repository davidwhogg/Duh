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

def populate_A_matrix(ones, image):
    """
    ## bugs:
    - Needs more information in this comment header.
    - Not yet written.
    - Ought to have user-settable order and options related to model complexity.
    - There must be a MUCH faster way to write this.
    """
    ny, nx = image.shape
    npix = nx * ny
    A = np.zeros((10, npix))
    A[0] = ones.reshape(npix)
    A[1] = image.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:, :-1] = image[:, 1:]
    A[2] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:, 1:] = image[:, :-1]
    A[3] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:-1] = image[1:]
    A[4] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[1:] = image[:-1]
    A[5] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:-1, :-1] = image[1:, 1:]
    A[6] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[:-1, 1:] = image[1:, :-1]
    A[7] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[1:, :-1] = image[:-1, 1:]
    A[8] = foo.reshape(npix)
    foo = np.zeros((ny, nx))
    foo[1:, 1:] = image[:-1, :-1]
    A[9] = foo.reshape(npix)
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
    A = np.vstack((populate_A_matrix(np.ones((ny, nx)), basis_image),
                   populate_A_matrix(xg, xg * basis_image),
                   populate_A_matrix(yg, yg * basis_image)))
    npix = ny * nx
    W = weight_image.reshape(npix)
    ATA = np.dot(W[None, :] * A, A.T)
    ATb = np.dot(W[None, :] * A, image_to_synthesize.reshape(npix))
    pars = np.linalg.solve(ATA, ATb)
    return np.dot(A.T, pars).reshape((ny, nx))

def synthesize_and_plot(fn1, fn2, fn):
    data_j2c = pf.open(fn1)[0].data
    data_j3c = pf.open(fn2)[0].data
    weight = np.ones_like(data_j2c)
    synth_j2c_j3c = synthesize_image(data_j2c, weight, data_j3c)
    resid_j2c_j3c = data_j2c - synth_j2c_j3c

    plt.clf()
    plt.subplot(221)
    plt.gray()
    vmin_j2c = np.percentile(data_j2c, 5.)
    vmax_j2c = np.percentile(data_j2c, 95.)
    y1 = 300
    y2 = y1 + 200
    x1 = 300
    x2 = x1 + 200
    plt.imshow(data_j2c[y1: y2, x1: x2], interpolation="nearest", vmin=vmin_j2c, vmax=vmax_j2c)
    plt.title(fn1)

    plt.subplot(222)
    vmin_j3c = np.percentile(data_j3c, 5.)
    vmax_j3c = np.percentile(data_j3c, 95.)
    plt.imshow(data_j3c[y1: y2, x1: x2], interpolation="nearest", vmin=vmin_j3c, vmax=vmax_j3c)
    plt.title(fn2)

    plt.subplot(223)
    plt.imshow(synth_j2c_j3c[y1: y2, x1: x2], interpolation="nearest", vmin=vmin_j2c, vmax=vmax_j2c)
    plt.title("synthetic " + fn1)

    plt.subplot(224)
    diff = np.median(data_j2c) - np.median(resid_j2c_j3c)
    plt.imshow(resid_j2c_j3c[y1: y2, x1: x2], interpolation="nearest", vmin=vmin_j2c - diff, vmax=vmax_j2c - diff)
    plt.title("residual")
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
