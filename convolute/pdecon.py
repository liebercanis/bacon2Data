import numpy as np
from scipy.signal import deconvolve

# Example: Signal with convolution
original_signal = np.array([1, 2, 3, 4, 5])
filter_kernel = np.array([1, 0.5])
convolved_signal = np.convolve(original_signal, filter_kernel)

# Deconvolution
recovered_signal, remainder = deconvolve(convolved_signal, filter_kernel)

print("Original Signal: ", original_signal)
print("Convolved Signal: ", convolved_signal)
print("Recovered Signal: ", recovered_signal)
