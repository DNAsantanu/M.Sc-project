# M.Sc-project
This repo contains some of the useful results and codes which I used during my M.Sc project (Project title: SKA Data challenge 3) in IIT Madras, Chennai

# 🛰️ SKA Data Challenge 3 – 21-cm Cosmological Signal Estimation

> A comprehensive project focused on estimating the 21-cm power spectrum using the Tapered Gridded Estimator (TGE), utilizing simulated interferometric data from the SKA1-Low telescope.

---

## 📘 Abstract

The 21-cm signal emitted by neutral hydrogen is a fundamental probe of the early universe, including the cosmic dawn and epoch of reionization. In this project, we validate the use of the **Tapered Gridded Estimator (TGE)** on data from the **SKA1-Low** telescope, aiming to recover both angular and cylindrical power spectra. Using interferometric radio data, we validate TGE against known models and explore implications for large-scale cosmological structure detection.

---

## 👨‍🔬 Contributor

- **Santanu Das** (PH22C040)  
  Department of Physics  
  Indian Institute of Technology Madras  
  **Guided by:** Dr. Samir Choudhuri

---

## 🌌 Scientific Motivation

The 21-cm line, caused by the hyperfine transition of neutral hydrogen, provides a unique probe into:

- The Dark Ages (z ≈ 200–25)
- Cosmic Dawn (z ≈ 25–16)
- Epoch of Reionization (z ≈ 10–6)
- Formation of first stars and galaxies
- Structure of primordial hydrogen clouds

Observing redshifted 21-cm radiation helps map these phases across cosmic time.

---

## 📡 Telescope & Data

### 🔭 SKA1-LOW Overview

- **Location**: Western Australia
- **Baseline**: Up to 65 km (131,072 antennas)
- **Sensitivity**: 400,000 m² effective area
- **Frequency Range**: 50–350 MHz (z=27–3)
- **Beamforming**: Digital; provides multiple simultaneous beams

### 📦 Dataset Description

- **Main Data**: 900 `uvfits` files (~7.5TB)
- **Test Data**: 150 `uvfits` files (166–181 MHz)
- **Station Beam Image File**: 22.4 GB (5°×5° FOV)
- **Channel Width**: 100 kHz
- **Observation Time**: 4 hours | **Integration Time**: 10 s

---

## 📖 Theoretical Background

### ✨ 21-cm Line Physics

- Arises from the spin-flip transition in neutral hydrogen.
- Observed frequency:  
   f_{observed} = 1420\(1 + z) MHz

### 📈 Visibility and Sky Image

- Sky brightness \( I(l,m) \) is Fourier-transformed to visibilities \( V(u,v) \) on the baseline plane.

### 📡 Primary Beam Pattern

- Directional sensitivity of each antenna affects the observed signal:
  \[
  V(u, v) = \int \int I(l,m) B(l,m) e^{-2\pi i (lu + mv)} dl dm
  \]

### 🌀 Tapered Gridded Estimator (TGE)

TGE mitigates sidelobe noise and reduces computation by convolving visibilities with a tapering function.  
Angular Power Spectrum \( C_\ell \) and Cylindrical Spectrum \( P(k_\perp, k_\parallel) \) are estimated from gridded visibilities.

---

## ⚙️ Data Processing Pipeline

### 🧩 Preprocessing Steps

- **CASA**: `split` to extract 1.5 km uv-region, `ms transform`, averaging
- Combined test data and realizations to create 10 FITS files
- TGE applied to:
  - **Test Data**: validation with known model
  - **True Data (Band 1: 106–121 MHz)**

### 📐 Estimation

- Angular Power Spectrum \( C_\ell(\Delta \nu) \)  
- Cylindrical Power Spectrum \( P(k_\perp, k_\parallel) \)

> Notably, test data show decorrelation with frequency—true data is noise dominated.

---

## 🔬 Key Results

| Feature | Observations |
|--------|-------------|
| **Validation** | Power spectrum matches well for low \( k_\perp, k_\parallel \) |
| **True Data** | Noise dominates power spectrum |
| **Efficiency** | Multi-script parallelization for faster averaging |
| **Visibility** | Phase/amplitude analyzed on uv-plane |

---

## 📊 Visual Outputs

- Amplitude & phase plots of gridded visibilities
- Angular power spectrum vs frequency separation
- Cylindrical power spectrum heatmaps
- Percentage deviation maps between model and output spectrum

---

## 📚 References

- Pritchard & Loeb, *21-cm Cosmology*, Reports on Progress in Physics
- Choudhuri et al., *TGE Estimator*, MNRAS (2016, 2018)
- Thompson & Moran, *Interferometry in Radio Astronomy*
- Braun et al., *SKA1 Anticipated Performance*

---

## 📬 Contact

**Santanu Das**  
Email: ph22c040@smail.iitm.ac.in  
M.Sc. Physics – IIT Madras

> 📁 This repository contains TGE implementation code and visualization scripts only. Raw data is confidential and not shared.

