import numpy as np
from numpy.fft import fft2, ifft2, fftshift, ifftshift
import matplotlib.pyplot as plt

# Create example image and kernel
image = np.zeros((100, 100))
image[40:60, 40:60] = 1
psf = np.zeros_like(image)
psf[48:52, 48:52] = 1

# Convolve using FFT
blurred = ifft2(fft2(image) * fft2(psf)).real

# Deconvolution
eps = 1e-6  # To avoid division by zero
deconvolved = ifft2(fft2(blurred) / (fft2(psf) + eps)).real

# Display results
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
axes[0].imshow(image, cmap='gray')
axes[0].set_title('Original Image')
axes[1].imshow(blurred, cmap='gray')
axes[1].set_title('Blurred Image')
axes[2].imshow(deconvolved, cmap='gray')
axes[2].set_title('Deconvolved Image')
plt.show()

