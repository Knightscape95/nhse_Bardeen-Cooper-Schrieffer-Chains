
# Non-Hermitian Quantum Metrology Enhancement
## PT-Symmetric Bardeen-Cooper-Schrieffer Chains

**Author:** Harshank Matkar  
**Institution:** Government College of Engineering Aurangabad  
**Email:** be22f03f033@geca.ac.in  
**Date:** August 24, 2025  

**Paper:** "Non-Hermitian Quantum Metrology Enhancement and Skin Effect Suppression in PT-Symmetric Bardeen-Cooper-Schrieffer Chains"

---

## Overview

This MATLAB codebase implements a comprehensive theoretical framework for quantum metrology in non-Hermitian systems, demonstrating:

1. **NHSE Suppression:** F_Q ∝ N³e^(-2κN) - exponential sensitivity failure
2. **PT Enhancement:** F_Q ∝ N²/δ - Heisenberg-limited enhancement  
3. **Multiparameter QFI:** Simultaneous estimation of (μ, φ, g)
4. **Experimental protocols** for superconducting circuits

## Quick Start Guide

### 1. Generate Data
```
main_data_generation()
% Creates comprehensive validation datasets
% Output: results_YYYYMMDD_HHMMSS/ folder with CORE_DATA_MASTER.mat
```

### 2. Generate Figures
```
run_complete_analysis_1()
% Produces all manuscript and supplementary figures
% Output: Enhanced analysis with publication-ready plots
% Present in config directory
```

### 3. Validate Functions
```
runtests('test_*.m')
% Run unit tests for all core functions
```

## Key Theoretical Results

**Enhancement Factors** (N=50, realistic parameters):
- **Chemical potential:** η_μ ≈ 20√N = 141
- **Peierls phase:** η_φ ≈ t²√(3N/2) = 100√N  
- **Gain/loss:** η_g ≈ Δ√N/(2t)

**QFI Matrix** (Heisenberg scaling F ∝ N²):
```
F = [N²/(4Δ²)    -N²Δt²/4      0    ]
    [-N²Δt²/4    3N²t⁴/2      0    ]  
    [0           0            N²Δ²/(4t²)]
```



## System Requirements
It was created on MATLAB Online.

**MATLAB Version:** R2020a or later

**Required Toolboxes:**
- Symbolic Math Toolbox (analytical calculations)
- Statistics and Machine Learning Toolbox (QFI analysis) 
- Parallel Computing Toolbox (optional, large-scale simulations)


## Workflow Overview

### 1. Parameter Configuration
- `load_params('pt')` → PT-symmetric systems
- `load_params('nhse')` → Skin effect models
- `load_params('experimental')` → Realistic circuits

### 2. Data Generation Pipeline
`main_data_generation()` creates:
- Exceptional points data (N = 4-50)
- Biorthogonal QFI systems (N = 4-100) 
- Multiparameter matrices (N = 10-200)
- PT-symmetric configurations (N = 6-200)
- NHSE suppression analysis (N = 8-50)

### 3. Analysis & Validation
`run_complete_analysis_1()` performs:
- NHSE exponential suppression verification
- PT exceptional point enhancement validation
- Multiparameter Heisenberg scaling confirmation
- Advanced theoretical predictions testing
- Publication figure generation

## Key Functions Usage Examples

### Example 1: Basic PT Enhancement Analysis
```
params = load_params('pt', 'N', 50);
results = pt_symmetric_bdg_analysis(params);
fprintf('Enhancement factor: η = %.1f\n', results.enhancement_factor);
```

### Example 2: Multiparameter QFI Matrix
```
[F, F_inv] = multiparameter_qfi(params);
disp('QFI Matrix (Heisenberg scaling N²):');
disp(F);
```

### Example 3: Exceptional Points Analysis  
```
ep_results = exceptional_points('all');
disp('PT-breaking thresholds for all systems');
```

### Example 4: NHSE Suppression Validation
```
nhse_results = nhse_analysis();
summarize_nhse_results(nhse_results);
```

## Experimental Parameters

**Superconducting Circuits** (from manuscript Sec. VI):
- **Hopping:** t/2π = 10 MHz (inductive coupling)
- **Pairing:** Δ/2π = 1 MHz (parametric drive)  
- **Gain/loss:** g/2π = 1 MHz (pump modulation)
- **System size:** N = 20-50 sites (modular design)
- **Decoherence:** γ/2π ∼ 50 kHz (T₂* = 100 μs)

**Expected Performance:**
- **Chemical potential sensing:** δμ ∼ 10⁻³/√ν
- **Phase estimation:** δφ ∼ 10⁻⁶/(Nt²√ν)
- **Enhancement over classical:** 10²-10³

## Validation & Benchmarking

**Numerical Accuracy Targets:**
- **Relative error:** < 10⁻¹² (analytical vs numerical)
- **Eigenvalue precision:** < 10⁻¹⁰ 
- **Biorthogonal overlap:** ||⟨ψ_L|ψ_R⟩ - δ_mn|| < 10⁻¹²
- **Matrix inversion:** ||F·F⁻¹ - I|| < 10⁻¹⁰

**Performance Benchmarks:**
- **Small systems** (N ≤ 20): Dense methods, ~seconds
- **Medium systems** (N = 20-100): Sparse methods, ~minutes  
- **Large systems** (N > 100): Parallel processing, ~hours
- **Memory scaling:** O(N²) dense, O(N) sparse

## Troubleshooting

### Common Issues & Solutions

1. **"gN ≥ 1 violates perturbative condition"**
   - **Solution:** Adjust g in params_pt.m to ensure gN < 0.95

2. **"Too close to exceptional point"**  
   - **Solution:** Increase safety_margin in PT configuration

3. **"Biorthogonal overlap too small"**
   - **Solution:** Check for NHSE regime, increase system size

4. **"Matrix inversion unstable"**
   - **Solution:** Use 'numerical' mode instead of 'analytical'

5. **"Memory insufficient for large N"**
   - **Solution:** Enable sparse methods in computational_settings.m

## Output Files & Data Structure

### main_data_generation() creates:
```
results_YYYYMMDD_HHMMSS/
├── CORE_DATA_MASTER.mat          # Complete dataset (primary)
├── exceptional_points_data.mat    # EP configurations
├── biorthogonal_qfi_data.mat     # QFI systems
├── multiparameter_qfi_data.mat   # Parameter matrices
├── pt_symmetric_bdg_data.mat     # BdG systems
├── nhse_analysis_data.mat        # NHSE suppression
├── VALIDATION_SUMMARY.txt        # Results summary
└── CORE_DATA_USAGE_GUIDE.txt     # Usage instructions
```

### run_complete_analysis_1() creates:
```
results_YYYYMMDD_HHMMSS/
├── COMPLETE_ANALYSIS_RESULTS.mat # Analysis results
├── ANALYSIS_SUMMARY.txt          # Performance report
└── figures/                      # Publication figures
    ├── Fig1_dichotomy.eps        # NHSE vs PT enhancement  
    ├── Fig2_multiparameter.eps   # Heisenberg scaling
    ├── FigS1_localization.eps    # Supplementary figures
    └── ...
```

## Citation

If you use this code in your research, please cite:

```
@article{matkar2025nonhermitian,
  title={Non-Hermitian Quantum Metrology Enhancement and Skin Effect 
         Suppression in PT-Symmetric Bardeen-Cooper-Schrieffer Chains},
  author={Matkar, Harshank},
  journal={https://arxiv.org/abs/2508.04815},
  year={2025},
  note={MATLAB code: DOI to be assigned}
}
```

## License

This project is licensed under the MIT License.  
Copyright (c) 2025 Harshank Matkar

## Support & Contributions

**For questions, bug reports, or contributions:**
- **Email:** be22f03f033@geca.ac.in
- **Issues:** Create detailed bug reports with system info
- **Contributions:** Fork repository and submit pull requests

Please include MATLAB version, system specifications, and error messages when reporting issues.

## Version History

**v1.0 (August 24, 2025):**
- Initial release with complete theoretical framework
- All manuscript results validated and reproducible
- Comprehensive test suite and documentation
- Ready for peer review and publication

---

## Getting Started

**To get started immediately:**

1. **Run:** `main_data_generation()`     → Generate validation datasets
2. **Run:** `run_complete_analysis_1()`  → Create figures and analysis  
3. **Check:** `results_*/ANALYSIS_SUMMARY.txt` for performance report

For detailed usage examples, see individual function documentation.
