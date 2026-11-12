import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, quad
from scipy.optimize import curve_fit
from scipy.stats import moyal, norm
from scipy.signal import convolve
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import os

# Create directory for saving plots
os.makedirs('BTB', exist_ok=True)

PLOT = True  
PLOT_Neu = True  
PLOT_ABS = True  

# =========================================================
# 1. Load Absorption Data & Perform Global Cross-Section Fit
# =========================================================
print("Loading data and fitting Gaussian cross-section to Neumeier data...")

c_ref, x_ref = 0.1, 11.6

try:
    data_rayleigh = np.genfromtxt('Neumeier_Rayleigh.csv', delimiter=',', skip_header=1)
    wl_rayleigh, t_rayleigh = data_rayleigh[:, 0], data_rayleigh[:, 1]
except OSError:
    print("Warning: Neumeier_Rayleigh.csv not found. Using dummy 0.90 flat baseline.")
    wl_rayleigh = np.linspace(110, 150, 100)
    t_rayleigh = np.full_like(wl_rayleigh, 0.90)

get_rayleigh = interp1d(wl_rayleigh, t_rayleigh, kind='linear', fill_value='extrapolate')

files = {
    0.1: 'Neumeier_absorption.csv',
    1.0: '1ppmabsorption.csv',
    10.0: '10ppmabsorption.csv',
    1000.0: '1000ppmabsorption.csv'
}

all_lam = []
all_c = []
all_trans_norm = []
data_dict = {}

C_DF = -0.025

for c_val, filename in files.items():
    try:
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    except OSError:
        wl_dummy = np.linspace(120, 132, 100)
        dummy_dip = 0.5 * np.exp(-((wl_dummy - 126.7)**2 / (2 * 1.0**2)))
        sim_trans = get_rayleigh(wl_dummy) * (np.clip(1.0 - dummy_dip, 1e-9, 1.0))**(c_val / 0.1)
        data = np.column_stack((wl_dummy, sim_trans + np.random.normal(0, 0.02, 100)))

    wl_raw = data[:, 0]
    trans_raw = data[:, 1]
    
    rayleigh_at_wl = get_rayleigh(wl_raw)
    trans_norm = trans_raw / (rayleigh_at_wl + C_DF)
#    trans_norm = trans_raw / (rayleigh_at_wl)
    
    mask = (wl_raw >= 118) & (wl_raw <= 132)
    wl_masked = wl_raw[mask]
    trans_norm_masked = trans_norm[mask]
    
    data_dict[c_val] = (wl_masked, trans_norm_masked)
    all_lam.extend(wl_masked)
    all_c.extend(np.full_like(wl_masked, c_val))
    all_trans_norm.extend(trans_norm_masked)

X_data_fit = np.vstack((np.array(all_lam), np.array(all_c)))
y_data_fit = np.array(all_trans_norm)

def convolved_gaussian_transmission(X, A1, E1, sigma1, A2, E2, sigma2, slit_sigma):
    lam_data = X[0]
    c_data = X[1]
    
    lam_hr = np.linspace(110, 140, 3000)
    dlam = lam_hr[1] - lam_hr[0]
    E_hr = 1239.84 / lam_hr 
    
    cs_1 = A1 * np.exp(-((E_hr - E1)**2) / (2 * sigma1**2))
    cs_2 = A2 * np.exp(-((E_hr - E2)**2) / (2 * sigma2**2))
    cs_total = cs_1 + cs_2
    
    slit_pixels = slit_sigma / dlam
    result = np.zeros_like(lam_data)
    
    for c_val in np.unique(c_data):
        T_true = np.exp(-cs_total * c_val * x_ref)
        T_conv = gaussian_filter1d(T_true, slit_pixels, mode='nearest')
        
        mask_c = (c_data == c_val)
        result[mask_c] = np.interp(lam_data[mask_c], lam_hr, T_conv)
        
    return result

initial_guesses = [50.0, 9.785, 0.05, 20.0, 10.16, 0.08, 0.6]
bounds = (
    [0.0, 9.60, 0.001, 0.0, 10.00, 0.001, 0.01], 
    [10000.0, 9.90, 0.500, 10000.0, 10.40, 0.500, 2.0]
)

popt, pcov = curve_fit(convolved_gaussian_transmission, X_data_fit, y_data_fit, p0=initial_guesses, bounds=bounds)

A1_cs, E1_cs, sigma1_cs, A2_cs, E2_cs, sigma2_cs, slit_fit = popt

print(f"Fit Complete. Slit Sigma: {slit_fit:.2f} nm")
print(f"A1_cs = {A1_cs:.2f}, E1_cs = {E1_cs:.2f}, sigma1_cs = {sigma1_cs:.2f}")
print(f"A2_cs = {A2_cs:.2f}, E2_cs = {E2_cs:.2f}, sigma2_cs = {sigma2_cs:.2f}")

# Temporarily increase absorption

A1_cs = A1_cs * 1.0
A2_cs = A2_cs * 1.0

if PLOT_ABS:
    plt.figure(figsize=(10, 6))
    colors = {0.1: 'black', 1.0: 'blue', 10.0: 'green', 1000.0: 'red'}
    lam_smooth = np.linspace(118, 132, 800)

    for c_val in files.keys():
        if c_val in data_dict:
            wl_norm, trans_norm = data_dict[c_val]
            plt.plot(wl_norm, trans_norm, '.', color=colors.get(c_val, 'black'), alpha=0.4, label=f'{c_val} ppm Data')
            
            c_array_smooth = np.full_like(lam_smooth, c_val)
            X_smooth = np.vstack((lam_smooth, c_array_smooth))
            trans_fit = convolved_gaussian_transmission(X_smooth, *popt)
            plt.plot(lam_smooth, trans_fit, '-', color=colors.get(c_val, 'black'), linewidth=2, label=f'{c_val} ppm Convolved Fit')

    plt.axhline(1.0, color='gray', linestyle='--', alpha=0.5, label='Baseline')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Normalized Transmission')
    plt.title('Slit-Convolved Gaussian Fit of Xenon Transmission')
    plt.xlim(118, 132)
    plt.ylim(-0.05, 1.05)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('BTB/convolved_gaussian_transmission_fit.jpg', dpi=300)
#    plt.show()
    plt.close()

# =========================================================
# 2. Initialization and Constants
# =========================================================
ppmx = np.array([0.0001, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0])
distance = np.array([1.0, 10.0, 20.0, 30.0])
ppmn = 0.0

taus = 4 # {FIX}
taut = 1184 # {FIX} based on undoped
taumix_T = 4700 # Should be 4700?  {FIX}
taumix_S = 35  # This is based on 20 - 50 ns and should probably {VARY} but {FIX} for undoped fits.
taumix_R = 250   # This is the relaxation time from singlet to triplet {VARY} but {FIX} for undoped fits. 
# The previous two time are estimated and should be fitted to the data at low doping levels
taux = 20 # {FIX}
tauquench = 2250 # We don't really know these numbers {VARY}
tauquenchmix = 3200 #  {VARY}? or just set equal to tauquench but {FIX} for undoped fits.

Dab = 2.2e-5  # I had 1.53e-5, but excimer should be between 2.0 and 2.5e-5 {FIX}
RArAr = 5.21e-8 # Half the bond length of Ar2 plus Ar VdW radius + Xe VdW radius {FIX}
RArArsinglet = 52.0e-8 # Foerster Resonant Energy Transfer radius, only happens for the singlet (dipole) state (just guess for now) {VARY}?
RArXe = 6.18e-8 # Should be 6.18e-8 {FIX}
RArN = 4.41e-8 # {FIX}
Avogadro = 6.022e23
rhoAr = 1.3954 #g/cm^3
molmassAr = 39.948 #g/mol
mrhoAr = rhoAr / molmassAr #mol/cm^3
NTotal = 36000

# Radiative cascade branching fraction:
# Fraction of 128 nm photons absorbed at 9.78 eV (n=2 Wannier-Mott / 1P1)
# that relax into resonant MS (3P1) vs metastable MT (3P2)
b_rad = 0.70  # {VARY} typically 0.50 - 0.85

# 3-Component Visible / NIR Neutral Bremsstrahlung & Recombination (per NTotal = 36,000 VUV)
# Total visible/NIR yield in pure LAr is ~200 - 500 ph/MeV (~0.5% - 1.2% of NTotal)
N_nbrs_prompt = NTotal * 0.0005  # ~20 ph/MeV: Pure fast NBrS (< 1 ns, tau_nbrs = 0.5 ns)
N_rec         = NTotal * 0.0020  # ~80 ph/MeV: Intermediate 36 ns recombination visible/NIR
f_vis_slow    = 0.0025           # ~100 ph/MeV: Slow visible/NIR component tracking T(t)

k_rec = 1.0 / 36.0               # {FIX} But {VARY} for undoped fits from SiPM 0
tau_nbrs = 0.5                   # Sub-nanosecond prompt NBrS thermalization time (ns)

R_singlet_triplet = 0.313  # {FIX} But {VARY} for undoped fits.
Sfrac = 1.0 / (1.0 + (1.0 / R_singlet_triplet))

t_min = 1300.0       # Pre-trigger baseline start (ns) to catch left Gaussian tail
t_start = 1450.0     # Physical excitation time t0 (ns)
tmax = 10000.0

d_target = 10.0
window_start = 1400.0  # Captures t_start - 5*sigma (1450 - 50 = 1400 ns)
window_end = 1510.0    # Captures t_start + 5*sigma + 3*taus (~1502 ns)

# Uniform 0.5 ns step ensures t_start = 1450.0 is an exact grid node
dt = 0.5
t_eval = np.arange(t_min, tmax + dt, dt)
post_mask = t_eval >= t_start
t_ode = t_eval[post_mask]  # Subset [1450.0, 10000.0] passed to the analytic solver

# =========================================================
# 3. ODE System Definition
# =========================================================
def lar_system(t, y, c_params):
    S, T, MS, MT, X, R = y
    (taus, taut, taudx_S, taudx_T, taudx_X, taudn, tauquench,
     taumix_S, taumix_T, taumix_R, tauquenchmix, taux, k_rec, Absorb, b_rad) = c_params
    
    # Total rate of 128 nm photons optically absorbed by Xe across the volume
    P_abs = Absorb * ((1.0 / taus) * S + (1.0 / taut) * T)
    
    dS  = -(1/taus + 1/taudx_S + 1/taudn + 1/tauquench) * S
    dT  = -(1/taut + 1/taudx_T + 1/taudn + 1/tauquench) * T
    dMS = b_rad * P_abs + (1/taudx_S)*S - (1/taumix_S + 1/taumix_R + 1/taudx_X + 1/taudn + 1/tauquenchmix) * MS
    dMT = (1.0 - b_rad) * P_abs + (1/taudx_T)*T + (1/taumix_R)*MS - (1/taumix_T + 1/taudx_X + 1/taudn + 1/tauquenchmix) * MT
    dX  = (1/taudx_X)*(MS + MT) - (1/taux + 1/taudn + 1/tauquench) * X
    dR  = -k_rec * R
    
    return [dS, dT, dMS, dMT, dX, dR]

# =========================================================
# 3b. Closed-Form (Analytic) Solution of lar_system
# =========================================================
# lar_system is a linear, constant-coefficient cascade: S and T each decay on
# their own (no incoming terms), MS is driven only by S and T, MT is driven
# only by S, T and MS, X is driven only by MS and MT, and R decays on its
# own. Because every equation's right-hand side depends only on states that
# are themselves already pure sums of exponentials, each state's solution is
# again a sum of exponentials, obtained via the standard integrating-factor
# identity for dy/dt + l*y = a*exp(-l'*t), y(0)=y0:
#   y(t) = y0*exp(-l*t) + a/(l-l') * (exp(-l'*t) - exp(-l*t))
# Repeating this substitution down the cascade (S,T -> MS -> MT -> X) gives
# every state as a "chain" of exponentials at the rates of every state
# upstream of it, plus its own rate:
#   chain(l1,...,ln; t) = sum_i exp(-li*t) / prod_{j!=i} (lj - li)
# which is exactly the n-fold time-convolution of exp(-l1*t),...,exp(-ln*t).
def chain(rates, t):
    """n-fold convolution of exp(-rates[i]*t); t may be a numpy array."""
    out = np.zeros_like(t, dtype=float)
    n = len(rates)
    for i in range(n):
        denom = 1.0
        for j in range(n):
            if j != i:
                denom *= (rates[j] - rates[i])
        out += np.exp(-rates[i] * t) / denom
    return out

def rates_from_params(c_params):
    """Build the lar_system rate constants from the same c_params tuple passed to lar_system."""
    (taus, taut, taudx_S, taudx_T, taudx_X, taudn, tauquench,
     taumix_S, taumix_T, taumix_R, tauquenchmix, taux, k_rec, Absorb, b_rad) = c_params

    lS = 1/taus + 1/taudx_S + 1/taudn + 1/tauquench
    lT = 1/taut + 1/taudx_T + 1/taudn + 1/tauquench
    lMS = 1/taumix_S + 1/taumix_R + 1/taudx_X + 1/taudn + 1/tauquenchmix
    lMT = 1/taumix_T + 1/taudx_X + 1/taudn + 1/tauquenchmix
    lX = 1/taux + 1/taudn + 1/tauquench

    # S and T each feed BOTH MS and MT, split by the radiative branching
    # fraction b_rad, via the shared absorbed-photon rate
    # P_abs = Absorb*((1/taus)*S + (1/taut)*T); S and T also each have a
    # separate direct diffusion feed (taudx_S -> MS, taudx_T -> MT) that is
    # NOT split by b_rad.
    aSS = b_rad * Absorb / taus + 1/taudx_S   # S -> MS
    aTS = b_rad * Absorb / taut                # T -> MS
    aST = (1 - b_rad) * Absorb / taus          # S -> MT
    aTT = (1 - b_rad) * Absorb / taut + 1/taudx_T  # T -> MT

    b = 1/taumix_R   # MS -> MT
    g = 1/taudx_X    # (MS, MT) -> X
    k = k_rec

    return (lS, lT, lMS, lMT, lX, aSS, aTS, aST, aTT, b, g, k)

def lar_analytic_solution(t, y0, rates):
    """Closed-form S,T,MS,MT,X,R at time(s) t (measured from the start of integration)."""
    S0, T0, MS0, MT0, X0, R0 = y0
    lS, lT, lMS, lMT, lX, aSS, aTS, aST, aTT, b, g, k = rates

    S = S0 * chain([lS], t)
    T = T0 * chain([lT], t)
    R = R0 * chain([k], t)

    MS = (MS0 * chain([lMS], t)
          + aSS * S0 * chain([lS, lMS], t)
          + aTS * T0 * chain([lT, lMS], t))

    MT = (MT0 * chain([lMT], t)
          + aST * S0 * chain([lS, lMT], t)
          + aTT * T0 * chain([lT, lMT], t)
          + b * MS0 * chain([lMS, lMT], t)
          + b * aSS * S0 * chain([lS, lMS, lMT], t)
          + b * aTS * T0 * chain([lT, lMS, lMT], t))

    X = (X0 * chain([lX], t)
         + g * (MS0 * chain([lMS, lX], t)
                + aSS * S0 * chain([lS, lMS, lX], t)
                + aTS * T0 * chain([lT, lMS, lX], t))
         + g * (MT0 * chain([lMT, lX], t)
                + aST * S0 * chain([lS, lMT, lX], t)
                + aTT * T0 * chain([lT, lMT, lX], t)
                + b * MS0 * chain([lMS, lMT, lX], t)
                + b * aSS * S0 * chain([lS, lMS, lMT, lX], t)
                + b * aTS * T0 * chain([lT, lMS, lMT, lX], t)))

    return S, T, MS, MT, X, R

# =========================================================
# 4. Convolution & Emission Profile Functions
# =========================================================
A_E, mu_E, sigma_E = 8132.0, 128.2, 3.006
Rayleigh_baseline = 1.00

def gaussian(x_val, amp, mu, sigma):
    return amp * np.exp(-((x_val - mu)**2 / (2*sigma**2))) 

def generate_emission_spectrum(wavelengths, y128, y150, y175, y500, x, c):
    mask = (wavelengths >= 110) & (wavelengths <= 900)
    shape_128 = np.zeros_like(wavelengths)
    lam_window = wavelengths[mask]
    
    em_128 = gaussian(lam_window, A_E, mu_E, sigma_E)
    
    E_window = 1239.84 / lam_window
    cs_1 = A1_cs * np.exp(-((E_window - E1_cs)**2) / (2 * sigma1_cs**2))
    cs_2 = A2_cs * np.exp(-((E_window - E2_cs)**2) / (2 * sigma2_cs**2))
    cs_total = cs_1 + cs_2
    
    t_scaled = np.exp(-cs_total * c * x) * Rayleigh_baseline
    shape_128[mask] = em_128 * t_scaled
    
    # Normalize 128nm shape so its Area = y128
    area_128 = trapezoid(shape_128, wavelengths)
    if area_128 > 0:
        shape_128 = (shape_128 / area_128) * y128

    # Widths of 150 and 175 nm excimer peaks
    sig_150 = 2.8
    sig_175 = 4.5
        
    amp_150 = y150 / (sig_150 * np.sqrt(2 * np.pi))
    amp_175 = y175 / (sig_175 * np.sqrt(2 * np.pi))

    # Broadband NBrS continuum (280 - 900 nm) + 760 nm Ar*(4p -> 4s) NIR hump
    nbrs_continuum = np.exp(-((wavelengths - 560.0) / 220.0)**4) * (wavelengths >= 250)
    nir_hump = 0.6 * gaussian(wavelengths, 1.0, 760.0, 55.0)
    shape_vis = nbrs_continuum + nir_hump
    area_vis = trapezoid(shape_vis, wavelengths)
    if area_vis > 0:
        shape_vis = (shape_vis / area_vis) * y500
        
    return (shape_128 +
            gaussian(wavelengths, amp_150, 149, sig_150) +
            gaussian(wavelengths, amp_175, 175, sig_175) +
            shape_vis)

def sipm_detector_qe(wavelengths):
    vuv4_wavelengths = np.array([
        120, 128, 150, 175, 200, 250, 300, 350, 400, 450, 
        500, 550, 600, 650, 700, 750, 800, 850, 900
    ])
    vuv4_pde = np.array([
        0.12, 0.15, 0.20, 0.24, 0.28, 0.33, 0.38, 0.44, 0.48, 0.50, 
        0.46, 0.38, 0.28, 0.20, 0.14, 0.09, 0.05, 0.03, 0.01
    ])
    get_vuv4_pde = interp1d(vuv4_wavelengths, vuv4_pde, kind='cubic', bounds_error=False, fill_value=0.0)
    return get_vuv4_pde(wavelengths)

def polycarbonate_trans(wavelengths):
    poly_wavelengths = np.array([
        120, 128, 150, 175, 200, 250, 300, 380, 410, 450, 
        500, 550, 600, 650, 700, 750, 800, 850, 900
    ])
    poly_trans = np.array([
        0.00, 0.00, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.85, 0.85, 
        0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85
    ])
    get_polytrans = interp1d(poly_wavelengths, poly_trans, kind='linear', bounds_error=False, fill_value=0.0)
    return np.clip(get_polytrans(wavelengths), 0.0, 1.0)

def quartz_trans(wavelengths):
    quartz_wavelengths = np.array([
        120, 130, 140, 150, 155, 160, 165, 170, 175, 180, 
        200, 300, 400, 500, 600, 700, 800, 850, 900
    ])
    quartz_trans = np.array([
        0.0036, 0.0022, 0.0014, 0.0035, .1344, .3840, .7008, .8605, .8862, .8946, 
        .9084, .9247, .9289, .9307, .9318, .9324, .9329, .9331, .9331
    ])
    get_quartztrans = interp1d(quartz_wavelengths, quartz_trans, kind='linear', bounds_error=False, fill_value=0.0)
    return np.clip(get_quartztrans(wavelengths), 0.0, 1.0)

def convolute_time_smear(t_eval, signal, sigma=10.0):
    dt = t_eval[1] - t_eval[0]
    half_pts = int(np.ceil(5.0 * sigma / dt))
    t_kernel = np.arange(-half_pts, half_pts + 1) * dt  # Odd length, exact 0.0 center
    kernel = norm.pdf(t_kernel, loc=0.0, scale=sigma)
    kernel /= np.sum(kernel)
    return convolve(signal, kernel, mode='same')

# =========================================================
# 5. Exact Integral Calculation
# =========================================================
def integrand_128(lam, x, c):
    em = gaussian(lam, A_E, mu_E, sigma_E)
    
    E = 1239.84 / lam
    cs_1 = A1_cs * np.exp(-((E - E1_cs)**2) / (2 * sigma1_cs**2))
    cs_2 = A2_cs * np.exp(-((E - E2_cs)**2) / (2 * sigma2_cs**2))
    cs_total = cs_1 + cs_2
    
    t_scaled = np.exp(-cs_total * c * x)
    return em * t_scaled

def calculate_total_light(x, c):
    total_light, _ = quad(integrand_128, 110, 140, args=(x, c))
    return total_light

I_0 = calculate_total_light(0.0, 0.0)

# =========================================================
# 6. Main Integration Loop
# =========================================================
# 2D arrays: rows = distance, columns = ppmx
window_integral_all = np.zeros((len(distance), len(ppmx)))
singlet_yield_all = np.zeros((len(distance), len(ppmx)))

for d_idx, davg in enumerate(distance):
    Total, Singlet, Triplet, Mixed, Xenon, Recomb = (np.zeros(len(ppmx)) for _ in range(6))
    Convoluted_Total = np.zeros(len(ppmx))

    for i, p in enumerate(ppmx):
        cX = 1e-6 * mrhoAr * Avogadro * p
        cN = 1e-6 * mrhoAr * Avogadro * ppmn
        
        taudx_T = 1.0 / (4 * np.pi * Dab * RArAr * cX * 1e-9)
        taudx_X = 1.0 / (4 * np.pi * Dab * RArXe * cX * 1e-9)
        taudx_S = 1.0 / (4 * np.pi * Dab * RArArsinglet * cX * 1e-9)
        taudn = 1.0 / (4 * np.pi * Dab * RArN * cN * 1e-9) if cN > 0 else np.inf
        
        Tr128 = calculate_total_light(davg, p) / I_0
        Absorb = 1.0 - Tr128
        
        y0 = [Sfrac * NTotal, (1 - Sfrac) * NTotal, 0, 0, 0, N_rec]
        params = (taus, taut, taudx_S, taudx_T, taudx_X, taudn, tauquench,
                  taumix_S, taumix_T, taumix_R, tauquenchmix, taux, k_rec, Absorb, b_rad)
        
        # Closed-form solution from t_start (1450 ns) to tmax (10000 ns);
        # t_ode - t_start is time since excitation, matching the analytic
        # solution's t'=0 convention
        rates = rates_from_params(params)
        S_a, T_a, MS_a, MT_a, X_a, R_a = lar_analytic_solution(t_ode - t_start, y0, rates)

        # Initialize full-grid arrays (zeros for t < 1450 ns)
        S_sol, T_sol, MS_sol, MT_sol, X_sol, R_sol = (np.zeros_like(t_eval) for _ in range(6))
        S_sol[post_mask], T_sol[post_mask], MS_sol[post_mask], MT_sol[post_mask], X_sol[post_mask], R_sol[post_mask] = S_a, T_a, MS_a, MT_a, X_a, R_a
        
        SLight = (1 - Absorb) * (1/taus) * S_sol
        TLight = (1 - Absorb) * (1/taut) * T_sol
        MLight = (1/taumix_S) * MS_sol + (1/taumix_T) * MT_sol
        XLight = (1/taux) * X_sol
        
        # 3-Component Visible/NIR Light (RLight):
        # 1. Prompt NBrS (< 1 ns) + 2. Intermediate Recomb (36 ns) + 3. Slow T-correlated visible/NIR
        prompt_nbrs = np.zeros_like(t_eval)
        prompt_nbrs[post_mask] = (N_nbrs_prompt / tau_nbrs) * np.exp(-(t_ode - t_start) / tau_nbrs)
        RLight = prompt_nbrs + k_rec * R_sol + f_vis_slow * (1/taut) * T_sol
        
        Singlet[i] = trapezoid(SLight, t_eval)
        Triplet[i] = trapezoid(TLight, t_eval)
        Mixed[i]   = trapezoid(MLight, t_eval)
        Xenon[i]   = trapezoid(XLight, t_eval)
        Recomb[i]  = trapezoid(RLight, t_eval)
        Total[i]   = Singlet[i] + Triplet[i] + Mixed[i] + Xenon[i] + Recomb[i]

        wavelengths = np.linspace(100, 900, 1000) 
        emission_curve = generate_emission_spectrum(wavelengths, (Singlet[i]+Triplet[i]), Mixed[i], Xenon[i], Recomb[i], davg, p)
        
        if PLOT:
            plt.figure(figsize=(8, 5))
            plt.plot(wavelengths, emission_curve, 'b', label='emission')
            plt.yscale('log')
            plt.ylim([0.005, 30000])
            plt.xlim([100, 600])
            plt.title(f"emission vs. wavelength, d = {davg}cm, ppmx = {p}")
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Yield')
            plt.legend()
            plt.savefig(f"BTB/emission_ppmx_{p}_davg_{davg}cm.jpg", dpi=300)
            plt.close()

        qe_curve = sipm_detector_qe(wavelengths)
        detected_spectrum = emission_curve * qe_curve
        Convoluted_Total[i] = trapezoid(detected_spectrum, wavelengths)

        mask128 = (wavelengths > 119) & (wavelengths < 140)
        SiPMQE128 = trapezoid(detected_spectrum[mask128], wavelengths[mask128])/trapezoid(emission_curve[mask128], wavelengths[mask128])
        mask150 = (wavelengths > 140) & (wavelengths < 160)
        SiPMQE150 = trapezoid(detected_spectrum[mask150], wavelengths[mask150])/trapezoid(emission_curve[mask150], wavelengths[mask150])
        mask175 = (wavelengths > 160) & (wavelengths < 180)
        SiPMQE175 = trapezoid(detected_spectrum[mask175], wavelengths[mask175])/trapezoid(emission_curve[mask175], wavelengths[mask175])
        mask500 = (wavelengths > 280) & (wavelengths <= 900)
        SiPMQE500 = trapezoid(detected_spectrum[mask500], wavelengths[mask500])/trapezoid(emission_curve[mask500], wavelengths[mask500])

        tot_time = (SiPMQE128*SLight + SiPMQE128*TLight + SiPMQE175*XLight + SiPMQE150*MLight + SiPMQE500*RLight)

        tot_time_conv = convolute_time_smear(t_eval, tot_time, sigma=10.0)
        slight_conv = convolute_time_smear(t_eval, SiPMQE128 * SLight, sigma=10.0)
        window_mask = (t_eval >= window_start) & (t_eval <= window_end)

        # Save into the row for the current distance (d_idx)
        window_integral_all[d_idx, i] = trapezoid(tot_time_conv[window_mask], t_eval[window_mask])
        singlet_yield_all[d_idx, i] = trapezoid(slight_conv[window_mask], t_eval[window_mask])

        if PLOT:
            plt.figure(figsize=(8, 5))
            plt.plot(t_eval, SiPMQE128 * SLight, 'b', label='Singlet')
            plt.plot(t_eval, SiPMQE128 * TLight, 'r', label='Triplet')
            plt.plot(t_eval, SiPMQE175 * XLight, 'g', label='Xenon')
            plt.plot(t_eval, SiPMQE150 * MLight, 'y', label='Mixed')
            plt.plot(t_eval, SiPMQE500 * RLight, 'm', label='Recomb')
            plt.plot(t_eval, tot_time, 'k', label='Total')
            plt.yscale('log')
            plt.ylim([0.005, 1000.0])
            duration = tmax - t_start
            plt.xlim([t_start - 0.1 * duration, t_start + 0.9 * duration])
            plt.title(f"SiPM yield vs. time, d = {davg}cm, ppmx = {p}")
            plt.xlabel('Time (ns)')
            plt.ylabel('Yield')
            plt.legend()
            plt.savefig(f"BTB/yield_time_SiPM_ppmx_{p}_davg_{davg}cm.jpg", dpi=300)
            plt.close()
            
            conv_signal = convolute_time_smear(t_eval, tot_time, 10)
            plt.figure(figsize=(8, 5))
            plt.plot(t_eval, conv_signal, 'k', label='Convoluted total')
            plt.yscale('log')
            plt.ylim([0.005, 1000.0])
            plt.xlim([t_start - 0.1 * duration, t_start + 0.9 * duration])
            plt.title(f"Convoluted SiPM yield vs. time, d = {davg}cm, ppmx = {p}")
            plt.xlabel('Time (ns)')
            plt.ylabel('Yield')
            plt.legend()
            plt.savefig(f"BTB/conv_yield_time_SiPM_ppmx_{p}_davg_{davg}cm.jpg", dpi=300)
            plt.close()

# --- NEW: Quartz-Filtered SiPM at d = 10.0 cm ---
            if np.isclose(davg, 10.0):
                q_curve = quartz_trans(wavelengths)
                det_spec_q = detected_spectrum * q_curve

                # Effective (SiPM QE * Quartz Transmission) across each emission band
                QE128_q = trapezoid(det_spec_q[mask128], wavelengths[mask128]) / trapezoid(emission_curve[mask128], wavelengths[mask128])
                QE150_q = trapezoid(det_spec_q[mask150], wavelengths[mask150]) / trapezoid(emission_curve[mask150], wavelengths[mask150])
                QE175_q = trapezoid(det_spec_q[mask175], wavelengths[mask175]) / trapezoid(emission_curve[mask175], wavelengths[mask175])
                QE500_q = trapezoid(det_spec_q[mask500], wavelengths[mask500]) / trapezoid(emission_curve[mask500], wavelengths[mask500])

                tot_time_q = (QE128_q * SLight + QE128_q * TLight +
                              QE150_q * MLight + QE175_q * XLight + QE500_q * RLight)
                conv_signal_q = convolute_time_smear(t_eval, tot_time_q, sigma=10.0)

                plt.figure(figsize=(8, 5))
                plt.plot(t_eval, conv_signal_q, 'g', linewidth=1.8, label='Convoluted total (Quartz)')
                plt.yscale('log')
                plt.ylim([1e-4, 1000.0])
                plt.xlim([t_start - 0.1 * duration, t_start + 0.9 * duration])
                plt.title(f"Convoluted Quartz-Filtered SiPM yield vs. time, d = {davg}cm, ppmx = {p}")
                plt.xlabel('Time (ns)')
                plt.ylabel('Yield')
                plt.legend()
                plt.savefig(f"BTB/conv_yield_time_Quartz_SiPM_ppmx_{p}_davg_{davg}cm.jpg", dpi=300)
                plt.close()

            # --- NEW: Polycarbonate-Filtered SiPM at d = 30.0 cm ---
            if np.isclose(davg, 30.0):
                p_curve = polycarbonate_trans(wavelengths)
                det_spec_p = detected_spectrum * p_curve

                # Effective (SiPM QE * Polycarbonate Transmission) across each emission band
                QE128_p = trapezoid(det_spec_p[mask128], wavelengths[mask128]) / trapezoid(emission_curve[mask128], wavelengths[mask128])
                QE150_p = trapezoid(det_spec_p[mask150], wavelengths[mask150]) / trapezoid(emission_curve[mask150], wavelengths[mask150])
                QE175_p = trapezoid(det_spec_p[mask175], wavelengths[mask175]) / trapezoid(emission_curve[mask175], wavelengths[mask175])
                QE500_p = trapezoid(det_spec_p[mask500], wavelengths[mask500]) / trapezoid(emission_curve[mask500], wavelengths[mask500])

                tot_time_p = (QE128_p * SLight + QE128_p * TLight +
                              QE150_p * MLight + QE175_p * XLight + QE500_p * RLight)
                conv_signal_p = convolute_time_smear(t_eval, tot_time_p, sigma=10.0)

                plt.figure(figsize=(8, 5))
                plt.plot(t_eval, conv_signal_p, 'm', linewidth=1.8, label='Convoluted total (Polycarbonate)')
                plt.yscale('log')
                # Lower ylim floor since Polycarbonate blocks VUV and only passes visible/NBrS
                plt.ylim([1e-5, 10.0])
                plt.xlim([t_start - 0.1 * duration, t_start + 0.9 * duration])
                plt.title(f"Convoluted Polycarbonate-Filtered SiPM yield vs. time, d = {davg}cm, ppmx = {p}")
                plt.xlabel('Time (ns)')
                plt.ylabel('Yield')
                plt.legend()
                plt.savefig(f"BTB/conv_yield_time_Polycarb_SiPM_ppmx_{p}_davg_{davg}cm.jpg", dpi=300)
                plt.close()

    if PLOT:
        plt.figure(figsize=(8, 5))
        plt.scatter(ppmx, Total, color='black', marker='p', label='Total (Static QE)')
        plt.scatter(ppmx, Convoluted_Total, color='purple', marker='*', s=100, label='Total (Convoluted QE)')
        plt.scatter(ppmx, Singlet, color='blue', marker='o', label='Singlet')
        plt.scatter(ppmx, Triplet, color='red', marker='s', label='Triplet')
        plt.scatter(ppmx, Mixed, color='gold', marker='d', label='Mixed')
        plt.scatter(ppmx, Xenon, color='green', marker='^', label='Xenon')
        plt.xscale('log')
        plt.xlim([0.00005, 50.0])
        plt.title(f"Light yield vs. Xenon conc., d = {davg}cm")
        plt.xlabel('Xenon Concentration (ppm)')
        plt.ylabel('Light Yield (NTotal = 30k)')
        plt.legend(loc='best')
        plt.savefig(f"BTB/yield_ppmx_davg_{davg}cm.jpg", dpi=300)
        plt.close()

# =========================================================
# 7. Post-Loop Plotting & Neumeier Replication
# =========================================================
if PLOT_Neu:
    t_eval = t_ode
    # Loop over every distance in your distance array
    for d_idx, davg in enumerate(distance):
        singlet_yield_current = singlet_yield_all[d_idx, :]
        window_integral_current = window_integral_all[d_idx, :]

        # Guard against division by zero if ppmx[0] is extremely small
        norm_factor_singlet = singlet_yield_current[0] if singlet_yield_current[0] != 0 else 1.0
        norm_factor_window = window_integral_current[0] if window_integral_current[0] != 0 else 1.0

        normalized_singlet = singlet_yield_current / norm_factor_singlet
        normalized_window = window_integral_current / norm_factor_window

        plt.figure(figsize=(8, 5))
        plt.plot(ppmx, normalized_singlet, 'bo-', linewidth=2, label='Normalized Total Singlet Yield (Model)')
        plt.plot(ppmx, normalized_window, 'r*-', linewidth=2, markersize=8, label=f'Normalized Convoluted Window ({window_start}-{window_end} ns)')
        plt.xscale('log')
        plt.xlabel('Xenon Concentration (ppm)')
        plt.ylabel('Normalized Light Yield')
        plt.title(f'Normalized Yield Comparison vs Xenon Concentration (d = {davg} cm)')
        plt.legend()
        plt.grid(True, which="both", ls="--", alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'BTB/normalized_windowed_yield_{davg}cm.jpg', dpi=300)
        plt.close()

    plt.figure(figsize=(8, 5))
    lam_window = np.linspace(110, 145, 500)
    em_128 = A_E * np.exp(-((lam_window - mu_E)**2) / (2 * sigma_E**2))
    normalization = 1.0 / max(em_128)
    em_128 = em_128 * normalization

    E_window = 1239.84 / lam_window
    cs_1 = A1_cs * np.exp(-((E_window - E1_cs)**2) / (2 * sigma1_cs**2))
    cs_2 = A2_cs * np.exp(-((E_window - E2_cs)**2) / (2 * sigma2_cs**2))
    cs_total = cs_1 + cs_2

    plt.plot(lam_window, em_128, 'k--', linewidth=2, label='Original 128 nm Emission')
    x_plot = d_target  # Uses d_target instead of hardcoded 10.0
    plot_concentrations = [0.0000001, 0.0001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0]

    for c_plot in plot_concentrations:
        t_scaled = np.exp(-cs_total * c_plot * x_plot) * Rayleigh_baseline
        surviving_128 = em_128 * t_scaled 
        plt.plot(lam_window, surviving_128, label=f'{c_plot} ppm Xenon')

    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Relative Intensity')
    plt.title(f'Surviving 128 nm Light Profile after Absorption (d = {x_plot} cm)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('BTB/surviving_128nm_shape.jpg', dpi=300)
    plt.close()

# --- NEW: Neumeier 2.3 mm Full Spectrum Replication Plot ---
    print("Generating Neumeier 2.3 mm Replication Plot...")
    plt.figure(figsize=(10, 6))
    
    neu_lam = np.linspace(115, 200, 2000)
    dlam_neu = neu_lam[1] - neu_lam[0]
    slit_pixels_neu = slit_fit / dlam_neu
    
    neu_concs = [0.0001, 0.1, 1.0, 10.0, 100.0, 1000.0]
    neu_colors = {0.0001: 'black', 0.1: 'black', 1.0: 'red', 10.0: 'green', 100.0: 'blue', 1000.0: 'magenta'}
    d_neu = 0.23 # 2.3 mm

    L_bulk = 4.0  # Average photon path length (cm) for the specific cell geometry

    for p in neu_concs:
        cX = 1e-6 * mrhoAr * Avogadro * p
        cN = 1e-6 * mrhoAr * Avogadro * ppmn
        
        taudx_T = 1.0 / (4 * np.pi * Dab * RArAr * cX * 1e-9)
        taudx_X = 1.0 / (4 * np.pi * Dab * RArXe * cX * 1e-9)
        taudx_S = 1.0 / (4 * np.pi * Dab * RArArsinglet * cX * 1e-9)
        taudn = 1.0 / (4 * np.pi * Dab * RArN * cN * 1e-9) if cN > 0 else np.inf

        # 1. Global Absorption (Feeds the ODEs to create M states)
        # Evaluated over the bulk dimension of the 2.3 mm cell
        Tr128_global = calculate_total_light(L_bulk, p) / I_0
        Absorb_global = 1.0 - Tr128_global

        # 2. Local Transmission (Feeds the final light yield plots)
        # Evaluated strictly over the path to the detector window
        Tr_detector = calculate_total_light(d_neu, p) / I_0
        
        y0 = [Sfrac * NTotal, (1 - Sfrac) * NTotal, 0, 0, 0, N_rec]
        params = (taus, taut, taudx_S, taudx_T, taudx_X, taudn, tauquench,
                  taumix_S, taumix_T, taumix_R, tauquenchmix, taux, k_rec, Absorb_global, b_rad)
        
        rates = rates_from_params(params)
        S_sol, T_sol, MS_sol, MT_sol, X_sol, R_sol = lar_analytic_solution(t_eval - t_start, y0, rates)
        
        SLight = Tr_detector * (1/taus) * S_sol
        TLight = Tr_detector * (1/taut) * T_sol
        MLight = (1/taumix_S) * MS_sol + (1/taumix_T) * MT_sol
        XLight = (1/taux) * X_sol
        
        prompt_nbrs_neu = (N_nbrs_prompt / tau_nbrs) * np.exp(-(t_eval - t_start) / tau_nbrs)
        RLight = prompt_nbrs_neu + k_rec * R_sol + f_vis_slow * (1/taut) * T_sol
        
        s_yield = trapezoid(SLight, t_eval)
        t_yield = trapezoid(TLight, t_eval)
        m_yield = trapezoid(MLight, t_eval) * 0.4 # Because much of the light is out of the acceptance of the spectrometer
        x_yield = trapezoid(XLight, t_eval)
        r_yield = trapezoid(RLight, t_eval)
        
        # Generates pure raw physical emission (No SiPM QE applied)
        raw_spectrum = generate_emission_spectrum(neu_lam, (s_yield+t_yield), m_yield, x_yield, r_yield, d_neu, p)
        
        # Convolve to match Neumeier's specific instrumental spectrometer blur
        convolved_spectrum = gaussian_filter1d(raw_spectrum, slit_pixels_neu, mode='nearest')
        
        plt.plot(neu_lam, convolved_spectrum, label=f'LAr {p} ppm Xe', color=neu_colors[p], linewidth=2)

    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Signal (Raw Simulated Photons)')
    plt.title('Simulated Raw Emission Spectrum at 2.3 mm (Neumeier Replication)')
    plt.xlim(115, 200)
    plt.legend(loc='upper right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('BTB/Neumeier_2.3mm_Replication.jpg', dpi=300)
    plt.show()
#    plt.close()
    print("Done!")
