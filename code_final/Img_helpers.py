import numpy as np
from scipy.ndimage import convolve

def bead_img(length, res, xs, ys, bead_sizes):
    num_pixels = int(length/res + 1)
    x = np.linspace(-length/2, length/2, num_pixels)
    y = np.linspace(-length/2, length/2, num_pixels)
    img = np.zeros((num_pixels, num_pixels))
    center = int(num_pixels/2)
    for k in range(len(xs)):
        x_center = xs[k]
        y_center = ys[k]
        bead_size = bead_sizes[k]
        for i in range(num_pixels):
            for j in range(num_pixels):
                x_img = (j - center) * (res)
                y_img = (i - center) * (res)

                distance = np.sqrt((x_center - x_img)**2 + (y_center - y_img)**2)
                if distance <= bead_size:
                    img[i, j] += 1
    #handle overlapping beads
    img[img > 1] = 1
    return x, y, img

def random_bead_img(length, res, num_beads, bead_sizes, seed = 1):
    rng = np.random.default_rng(seed)
    x_centers = rng.uniform(-length/2, length/2, num_beads)
    y_centers = rng.uniform(-length/2, length/2, num_beads)
    return bead_img(length, res, x_centers, y_centers, bead_sizes)

def tess(length, res, num_centers, seed = 1):
    rng = np.random.default_rng(seed)
    num_pixels = int(length/res + 1)
    x = np.linspace(-length/2, length/2, num_pixels)
    y = np.linspace(-length/2, length/2, num_pixels)
    x_centers = rng.uniform(-length/2, length/2, num_centers)
    y_centers = rng.uniform(-length/2, length/2, num_centers)
    img = np.zeros((num_pixels, num_pixels))
    center = int(num_pixels/2)
    for i in range(num_pixels):
        for j in range(num_pixels):
            best_center = -1 
            best_distance = np.inf
            for k in range(len(x_centers)):
                x_center = x_centers[k]
                y_center = y_centers[k]
                x_img = (j - center) * (res)
                y_img = (i - center) * (res)
                distance = np.sqrt((x_center - x_img)**2 + (y_center - y_img)**2)
                if distance < best_distance: 
                    best_center = k
                    best_distance = distance
            img[i ,j] = best_center
    return x, y, img

def voronoi(length, res, num_centers, seed =1):
    tessalation = 0
    x, y, tessalation  = tess(length, res, num_centers, seed)
    img = np.zeros_like(tessalation)
    num_pixels = int(length/res + 1)
    for i in range(1, num_pixels - 1):
        for j in range(1, num_pixels - 1):
            roi = tessalation[i-1:i+2, j-1:j+2]
            if len(np.unique(roi)) > 1:
                img[i, j] = 1
    return x, y, img

def cell(length, res, num_centers, seed = 1, pad_amt = 0):
    x, y, vor_img = voronoi(length, res, num_centers, seed)
    num_beads = num_centers * 5
    rng = np.random.default_rng(seed)
    _, _, bead_img = random_bead_img(length, res, num_beads, np.abs(rng.normal(loc = 0.0001, scale = 0.00005, size = num_beads)), seed)
    final_img = vor_img + bead_img 
    final_img[final_img > 1] = 1
    if pad_amt !=0:
        final_img[0:pad_amt, :] = 0
        final_img[-pad_amt:-1, :] = 0
        final_img[:, 0:pad_amt] = 0
        final_img[:, -pad_amt:-1] = 0
    return x, y, final_img

def psf_convolve(image, psf):
    return convolve(image, psf/np.sum(psf), mode="reflect")

def quality(image):
    return np.mean(image**2)
