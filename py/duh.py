"""
This file is part of the Duh project.
Copyright 2014 David W. Hogg (NYU) and Dustin Lang (CMU).

## bugs
* not yet written
* does nothing
"""

import numpy as np
import astropy.io.fits as pf

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

if __name__ == "__main__":
    new = np.random.normal(size=(1700, 1500))
    weight = np.ones((1700, 1500))
    old = np.random.uniform(size=(1700, 1500))
    print new - synthesize_image(new, weight, old)
