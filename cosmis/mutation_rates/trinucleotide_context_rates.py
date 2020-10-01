
MUTATION_RATES = {
     'AAA': {'ACA': 8.691e-06, 'AGA': 1.743e-05, 'ATA': 4.084e-06},
     'AAC': {'ACC': 7.53e-06, 'AGC': 3.313e-05, 'ATC': 7.305e-06},
     'AAG': {'ACG': 9.637e-06, 'AGG': 2.822e-05, 'ATG': 5.086e-06},
     'AAT': {'ACT': 6.682e-06, 'AGT': 4.623e-05, 'ATT': 6.913e-06},
     'ACA': {'AAA': 1.269e-05, 'AGA': 1.458e-05, 'ATA': 4.907e-05},
     'ACC': {'AAC': 3.334e-05, 'AGC': 2.001e-05, 'ATC': 6.951e-05},
     'ACG': {'AAG': 4.789e-05, 'AGG': 4.367e-05, 'ATG': 0.0007943},
     'ACT': {'AAT': 9.85e-06, 'AGT': 1.596e-05, 'ATT': 4.452e-05},
     'AGA': {'AAA': 3.053e-05, 'ACA': 2.124e-05, 'ATA': 1.089e-05},
     'AGC': {'AAC': 5.492e-05, 'ACC': 1.73e-05, 'ATC': 1.261e-05},
     'AGG': {'AAG': 6.116e-05, 'ACG': 2.095e-05, 'ATG': 1.057e-05},
     'AGT': {'AAT': 4.778e-05, 'ACT': 1.716e-05, 'ATT': 1.091e-05},
     'ATA': {'AAA': 7.792e-06, 'ACA': 5.365e-05, 'AGA': 6.003e-06},
     'ATC': {'AAC': 1.42e-05, 'ACC': 3.319e-05, 'AGC': 7.57e-06},
     'ATG': {'AAG': 1.221e-05, 'ACG': 8.107e-05, 'AGG': 1.215e-05},
     'ATT': {'AAT': 6.643e-06, 'ACT': 4.509e-05, 'AGT': 6.266e-06},
     'CAA': {'CCA': 1.216e-05, 'CGA': 2.756e-05, 'CTA': 5.441e-06},
     'CAC': {'CCC': 1.434e-05, 'CGC': 4.405e-05, 'CTC': 9.51e-06},
     'CAG': {'CCG': 1.348e-05, 'CGG': 5.227e-05, 'CTG': 9.337e-06},
     'CAT': {'CCT': 1.161e-05, 'CGT': 8.011e-05, 'CTT': 1.209e-05},
     'CCA': {'CAA': 1.2e-05, 'CGA': 1.368e-05, 'CTA': 4.878e-05},
     'CCC': {'CAC': 2.281e-05, 'CGC': 2.678e-05, 'CTC': 7.393e-05},
     'CCG': {'CAG': 7.323e-05, 'CGG': 0.0001038, 'CTG': 0.0009883},
     'CCT': {'CAT': 1.009e-05, 'CGT': 2.053e-05, 'CTT': 6.478e-05},
     'CGA': {'CAA': 0.0006176, 'CCA': 5.302e-05, 'CTA': 4.017e-05},
     'CGC': {'CAC': 0.0008322, 'CCC': 5.78e-05, 'CTC': 0.0001083},
     'CGG': {'CAG': 0.0009205, 'CCG': 9.615e-05, 'CTG': 6.193e-05},
     'CGT': {'CAT': 0.0007958, 'CCT': 4.554e-05, 'CTT': 4.364e-05},
     'CTA': {'CAA': 4.656e-06, 'CCA': 2.984e-05, 'CGA': 6.465e-06},
     'CTC': {'CAC': 6.6e-06, 'CCC': 2.303e-05, 'CGC': 1.017e-05},
     'CTG': {'CAG': 8.62e-06, 'CCG': 5.482e-05, 'CGG': 1.337e-05},
     'CTT': {'CAT': 5.364e-06, 'CCT': 2.902e-05, 'CGT': 9.346e-06},
     'GAA': {'GCA': 7.521e-06, 'GGA': 1.86e-05, 'GTA': 4.619e-06},
     'GAC': {'GCC': 8.077e-06, 'GGC': 3.265e-05, 'GTC': 1.384e-05},
     'GAG': {'GCG': 9.361e-06, 'GGG': 2.195e-05, 'GTG': 6.415e-06},
     'GAT': {'GCT': 7.022e-06, 'GGT': 3.081e-05, 'GTT': 1.304e-05},
     'GCA': {'GAA': 1.77e-05, 'GGA': 1.219e-05, 'GTA': 4.503e-05},
     'GCC': {'GAC': 2.603e-05, 'GGC': 2.137e-05, 'GTC': 7.488e-05},
     'GCG': {'GAG': 9.186e-05, 'GGG': 5.588e-05, 'GTG': 0.0007982},
     'GCT': {'GAT': 1.422e-05, 'GGT': 1.812e-05, 'GTT': 6.029e-05},
     'GGA': {'GAA': 4.937e-05, 'GCA': 1.959e-05, 'GTA': 1.387e-05},
     'GGC': {'GAC': 7.74e-05, 'GCC': 2.177e-05, 'GTC': 2.808e-05},
     'GGG': {'GAG': 7.727e-05, 'GCG': 2.823e-05, 'GTG': 2.164e-05},
     'GGT': {'GAT': 6.844e-05, 'GCT': 1.917e-05, 'GTT': 3.069e-05},
     'GTA': {'GAA': 6.772e-06, 'GCA': 3.752e-05, 'GGA': 6.628e-06},
     'GTC': {'GAC': 1.266e-05, 'GCC': 3.315e-05, 'GGC': 9.482e-06},
     'GTG': {'GAG': 9.078e-06, 'GCG': 4.386e-05, 'GGG': 1.268e-05},
     'GTT': {'GAT': 8.112e-06, 'GCT': 3.533e-05, 'GGT': 8.324e-06},
     'TAA': {'TCA': 7.661e-06, 'TGA': 2.133e-05, 'TTA': 9.146e-06},
     'TAC': {'TCC': 6.65e-06, 'TGC': 3.917e-05, 'TTC': 6.823e-06},
     'TAG': {'TCG': 6.676e-06, 'TGG': 2.899e-05, 'TTG': 5.558e-06},
     'TAT': {'TCT': 5.198e-06, 'TGT': 5.237e-05, 'TTT': 8.189e-06},
     'TCA': {'TAA': 8.46e-06, 'TGA': 1.389e-05, 'TTA': 3.162e-05},
     'TCC': {'TAC': 1.624e-05, 'TGC': 2.053e-05, 'TTC': 5.27e-05},
     'TCG': {'TAG': 3.495e-05, 'TGG': 6.161e-05, 'TTG': 0.0005842},
     'TCT': {'TAT': 1.201e-05, 'TGT': 2.194e-05, 'TTT': 3.078e-05},
     'TGA': {'TAA': 3.197e-05, 'TCA': 1.335e-05, 'TTA': 7.921e-06},
     'TGC': {'TAC': 4.539e-05, 'TCC': 1.289e-05, 'TTC': 1.761e-05},
     'TGG': {'TAG': 4.466e-05, 'TCG': 1.228e-05, 'TTG': 1.13e-05},
     'TGT': {'TAT': 5.093e-05, 'TCT': 1.367e-05, 'TTT': 1.311e-05},
     'TTA': {'TAA': 9.494e-06, 'TCA': 2.284e-05, 'TGA': 7.072e-06},
     'TTC': {'TAC': 4.846e-06, 'TCC': 1.973e-05, 'TGC': 7.974e-06},
     'TTG': {'TAG': 4.262e-06, 'TCG': 2.737e-05, 'TGG': 1.124e-05},
     'TTT': {'TAT': 3.996e-06, 'TCT': 1.737e-05, 'TGT': 8.759e-06}
}

MUTATION_RATES_UNIQUE = {
     'AAA': {'ACA': 4.38e-09, 'AGA': 7.47e-09, 'ATA': 2.51e-09},
     'AAC': {'ACC': 3.52e-09, 'AGC': 1.39e-08, 'ATC': 3.65e-09},
     'AAG': {'ACG': 4.1e-09, 'AGG': 1.05e-08, 'ATG': 2.55e-09},
     'AAT': {'ACT': 3.41e-09, 'AGT': 2.06e-08, 'ATT': 4.09e-09},
     'ACA': {'AAA': 7.44e-09, 'AGA': 5.81e-09, 'ATA': 2.34e-08},
     'ACC': {'AAC': 1.57e-08, 'AGC': 6.38e-09, 'ATC': 2.68e-08},
     'ACG': {'AAG': 1.13e-08, 'AGG': 7.55e-09, 'ATG': 2.35e-07},
     'ACT': {'AAT': 5.74e-09, 'AGT': 7.57e-09, 'ATT': 2.11e-08},
     'AGA': {'AAA': 1.42e-08, 'ACA': 8.58e-09, 'ATA': 7.17e-09},
     'AGC': {'AAC': 2.19e-08, 'ACC': 5.63e-09, 'ATC': 5.74e-09},
     'AGG': {'AAG': 2.32e-08, 'ACG': 6.55e-09, 'ATG': 5.01e-09},
     'AGT': {'AAT': 2.12e-08, 'ACT': 7.67e-09, 'ATT': 5.77e-09},
     'ATA': {'AAA': 4.84e-09, 'ACA': 2.84e-08, 'AGA': 3.13e-09},
     'ATC': {'AAC': 8.22e-09, 'ACC': 1.4e-08, 'AGC': 2.9e-09},
     'ATG': {'AAG': 6.63e-09, 'ACG': 3.12e-08, 'AGG': 5.32e-09},
     'ATT': {'AAT': 4.02e-09, 'ACT': 2.04e-08, 'AGT': 3.35e-09},
     'CAA': {'CCA': 5.06e-09, 'CGA': 1.07e-08, 'CTA': 2.58e-09},
     'CAC': {'CCC': 4.69e-09, 'CGC': 1.38e-08, 'CTC': 3.81e-09},
     'CAG': {'CCG': 4.08e-09, 'CGG': 1.41e-08, 'CTG': 3.15e-09},
     'CAT': {'CCT': 5.22e-09, 'CGT': 3.05e-08, 'CTT': 6.53e-09},
     'CCA': {'CAA': 6.18e-09, 'CGA': 4.39e-09, 'CTA': 1.84e-08},
     'CCC': {'CAC': 7.75e-09, 'CGC': 7.52e-09, 'CTC': 2.39e-08},
     'CCG': {'CAG': 1.01e-08, 'CGG': 9.68e-09, 'CTG': 2.25e-07},
     'CCT': {'CAT': 5e-09, 'CGT': 6.47e-09, 'CTT': 2.31e-08},
     'CGA': {'CAA': 1.91e-07, 'CCA': 1.03e-08, 'CTA': 1.03e-08},
     'CGC': {'CAC': 2.02e-07, 'CCC': 6.79e-09, 'CTC': 1.63e-08},
     'CGG': {'CAG': 2.11e-07, 'CCG': 9.09e-09, 'CTG': 9.44e-09},
     'CGT': {'CAT': 2.32e-07, 'CCT': 7.49e-09, 'CTT': 1.11e-08},
     'CTA': {'CAA': 2.92e-09, 'CCA': 1.39e-08, 'CGA': 3.02e-09},
     'CTC': {'CAC': 2.88e-09, 'CCC': 8.32e-09, 'CGC': 3.13e-09},
     'CTG': {'CAG': 3.27e-09, 'CCG': 1.47e-08, 'CGG': 4.25e-09},
     'CTT': {'CAT': 2.44e-09, 'CCT': 1.01e-08, 'CGT': 3.95e-09},
     'GAA': {'GCA': 2.94e-09, 'GGA': 7.65e-09, 'GTA': 2.36e-09},
     'GAC': {'GCC': 2.56e-09, 'GGC': 1.2e-08, 'GTC': 5.2e-09},
     'GAG': {'GCG': 3.13e-09, 'GGG': 8.22e-09, 'GTG': 2.85e-09},
     'GAT': {'GCT': 2.8e-09, 'GGT': 1.35e-08, 'GTT': 7.94e-09},
     'GCA': {'GAA': 9.45e-09, 'GGA': 4.12e-09, 'GTA': 1.76e-08},
     'GCC': {'GAC': 1.07e-08, 'GGC': 5.97e-09, 'GTC': 2.46e-08},
     'GCG': {'GAG': 1.61e-08, 'GGG': 6.7e-09, 'GTG': 1.97e-07},
     'GCT': {'GAT': 5.91e-09, 'GGT': 5.75e-09, 'GTT': 2.27e-08},
     'GGA': {'GAA': 1.94e-08, 'GCA': 6.08e-09, 'GTA': 6.45e-09},
     'GGC': {'GAC': 2.47e-08, 'GCC': 6.04e-09, 'GTC': 1.06e-08},
     'GGG': {'GAG': 2.37e-08, 'GCG': 7.49e-09, 'GTG': 7.68e-09},
     'GGT': {'GAT': 2.68e-08, 'GCT': 6.37e-09, 'GTT': 1.56e-08},
     'GTA': {'GAA': 3.2e-09, 'GCA': 1.73e-08, 'GGA': 2.66e-09},
     'GTC': {'GAC': 5.37e-09, 'GCC': 1.23e-08, 'GGC': 2.67e-09},
     'GTG': {'GAG': 3.71e-09, 'GCG': 1.35e-08, 'GGG': 4.59e-09},
     'GTT': {'GAT': 3.74e-09, 'GCT': 1.46e-08, 'GGT': 3.6e-09},
     'TAA': {'TCA': 3.84e-09, 'TGA': 1.15e-08, 'TTA': 4.5e-09},
     'TAC': {'TCC': 2.71e-09, 'TGC': 1.75e-08, 'TTC': 3.24e-09},
     'TAG': {'TCG': 2.98e-09, 'TGG': 1.36e-08, 'TTG': 2.87e-09},
     'TAT': {'TCT': 3.17e-09, 'TGT': 2.83e-08, 'TTT': 4.8e-09},
     'TCA': {'TAA': 5.37e-09, 'TGA': 5.19e-09, 'TTA': 1.49e-08},
     'TCC': {'TAC': 6.89e-09, 'TGC': 6.43e-09, 'TTC': 2.05e-08},
     'TCG': {'TAG': 9.69e-09, 'TGG': 9.49e-09, 'TTG': 1.75e-07},
     'TCT': {'TAT': 7.19e-09, 'TGT': 8.6e-09, 'TTT': 1.42e-08},
     'TGA': {'TAA': 1.5e-08, 'TCA': 5.21e-09, 'TTA': 5.31e-09},
     'TGC': {'TAC': 1.81e-08, 'TCC': 4.19e-09, 'TTC': 9.56e-09},
     'TGG': {'TAG': 1.76e-08, 'TCG': 4.2e-09, 'TTG': 5.91e-09},
     'TGT': {'TAT': 2.32e-08, 'TCT': 5.8e-09, 'TTT': 7.34e-09},
     'TTA': {'TAA': 4.53e-09, 'TCA': 1.16e-08, 'TGA': 3.86e-09},
     'TTC': {'TAC': 2.4e-09, 'TCC': 7.78e-09, 'TGC': 2.95e-09},
     'TTG': {'TAG': 2.46e-09, 'TCG': 1.02e-08, 'TGG': 4.76e-09},
     'TTT': {'TAT': 2.48e-09, 'TCT': 7.35e-09, 'TGT': 4.27e-09}
}