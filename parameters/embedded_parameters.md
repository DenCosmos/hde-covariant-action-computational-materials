# Live configurations and embedded parameters

`figure01.json`, `figure02.json` and `noncom_windows.json` are read by their drivers. In `noncom.json`, momenta and release metadata are read; the speeds, signs and norms document the unchanged engine and are checked against it rather than replacing its constants. Do not assume every documentary field changes a calculation.

The radial generator embeds symbolic r and t=tan(theta/2), U=0, g_R=k=1, scalar speed 1/sqrt(3), tensor speed 1 and signed momenta. The tensor normalization conventions are in REPRODUCIBILITY.md. The angular driver reads `--j` (default 0,2,3,4); the stored covariance is parameterized, not chosen numerically. There is no random number generator.

The inherited background diagnostic script does **not** read an external parameter JSON. It uses DOP853 and the following explicit constants in the named functions. Units are t H0, L H0, H/H0, k/H0, chi/(M_Pl² H0²); setting dimensionless values to one does not identify H0 with M_Pl.

`saddle_numerics`: n=10,12,14; parallel/transverse eigenvector directions; amplitudes 1e-4,1e-5,1e-6; final dimensionless time 0.5; rtol 1e-9 and 1e-12; atol=rtol*amplitude*0.01. The n=10,12,14 finest-amplitude maximum error must be below 1e-7.

`crossing_numerics`: n=12, Lc=1.02, ac=1, k=1, chi_initial=-1e-2; the zero-chi background fixes Hc from G(Lc). Integrate in log(-chi) down to 1e-12, 1501 samples; two fundamental perturbation initial conditions (ell,ell_dot)=(1,0),(0,1); rtol 1e-9 and 1e-12, atol=rtol*0.01. Comparison scales are 1e-4,1e-8,1e-12. The final logarithmic-slope/flux difference must be below 1e-8.

`late_numerics`: n=10 and 12, d0=0.1, k=1; final log(a) 8 or 5; two fundamental real solutions (1,0),(0,1); 1200 samples per unit log(a); rtol 1e-9 and 1e-12, atol=rtol*0.01. The n=10 slope windows are [4,5],[6,7],[7,8]; n=12 uses [2,3],[3,4],[4,5]. The manuscript transcription tolerance is 0.51 of the last displayed decimal unit, not an arbitrary percent.

Figure 3 embeds its original c sampling and normalized cubic estimates in `article_fig03_plot_strong_coupling.py`. It does not read these documentary notes. Original parameters are retained; the script writes the actual plotted values to CSV.

Exact symbolic tests use generic variables with assumptions stated next to their declarations. A finite sampled substitution is never labeled an all-variable identity. Consult the source and per-test map for the exact inequalities and denominator exclusions.
