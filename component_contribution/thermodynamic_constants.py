R = 8.31e-3 # kJ/(K*mol)
F = 96.485 # kC/mol
J_per_cal = 4.184
default_T = 298.15 # K
default_I = 0.25 # M
default_pH = 7.0
default_c0 = 1 # M
default_pMg = 10
default_RT = R * default_T
default_c_mid = 1e-3 # M
default_c_range = (1e-6, 1e-2) # M
dG0_f_Mg = -455.3 # kJ/mol, formation energy of Mg2+

symbol_d_G = "&Delta;G"
symbol_d_G0 = "&Delta;G&deg;"
symbol_d_G_prime = "&Delta;G'"
symbol_d_G0_prime = "&Delta;G'&deg;"

symbol_dr_G = "&Delta;<sub>r</sub>G"
symbol_dr_G0 = "&Delta;<sub>r</sub>G&deg;"
symbol_dr_G_prime = "&Delta;<sub>r</sub>G'"
symbol_dr_G0_prime = "&Delta;<sub>r</sub>G'&deg;"
symbol_dr_Gc_prime = "&Delta;<sub>r</sub>G'<sup>c</sup>"

symbol_df_G = "&Delta;<sub>f</sub>G"
symbol_df_G0 = "&Delta;<sub>f</sub>G&deg;"
symbol_df_G_prime = "&Delta;<sub>f</sub>G'"
symbol_df_G0_prime = "&Delta;<sub>f</sub>G'&deg;"

# Approximation of the temperature dependency of ionic strength effects
DH_alpha = lambda T : 1e-3*(9.20483*T) - 1e-5*(1.284668 * T**2) + 1e-8*(4.95199 * T**3)
DH_beta = 1.6

# Debye-Huckel
debye_huckel = lambda I, T : DH_alpha(T) * I**(0.5) / (1.0 + DH_beta * I**(0.5))

