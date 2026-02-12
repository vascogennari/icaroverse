# Available PSD for simulation pipeline

We list here the PSD files that are available in the present simulation pipeline for population & cosmology studies built around `icarogw` hierarchical inference python package. These PSDs are used for SNR computation (either optimal SNR with `pycbc` or matched filter SNR with `bilby`) and realistic single events parameter estimation of simulated population. 

The available files are taken from [Public LIGO DCC](https://dcc.ligo.org/LIGO-T2000012/public). NB: all the files actually contain ASD data (reminder : $\rm ASD = \sqrt{PSD}$). This is good to keep in mind, but the files are handled accordingly in the pipeline anyway.

### Available detectors & corresponding available observing run sensitivities

| Detector        | id  | O3  | O4  | O5  |
| :-------------- | :-: | :-: | :-: | :-: |
| LIGO Livingston | L1  | YES | YES | YES |
| LIGO Handford   | H1  | YES | YES | YES |
| Virgo           | V1  | YES | YES | YES |
| KAGRA           | K1  | NO  | YES | YES |

For each entry of this table, one or several files are provided. We list them below, with specific information.

| File name               | detector | run | BNS range          | comment                                  |
| :--------------         | :------- | :-- | :--------------    | :--------------------------------------- |
| `aligo_O3actual_L1.txt` | L1       | O3  | $130 ~ \rm Mpc$    | based on first three months of O3 |
| `aligo_O3actual_H1.txt` | H1       | O3  | $110 ~ \rm Mpc$    | based on first three months of O3 |
| `avirgo_O3actual.txt`   | V1       | O3  | $50 ~ \rm Mpc$     | based on first three months of O3 |
| `aligo_O4low.txt`       | H1, L1   | O4  | $160 ~ \rm Mpc$    | for O4 low sensitivity simulations |
| `aligo_O4high.txt`      | H1, L1   | O4  | $190 ~ \rm Mpc$    | for O4 high sensitivity simulations |
| `avirgo_O4high_NEW.txt` | V1       | O4  | $90-120 ~ \rm Mpc$ | for O4 simulations |
| `kagra_3Mpc.txt`        | K1       | O4  | $3 ~ \rm Mpc$      | for O4 low sensitivity simulations |
| `kagra_10Mpc.txt`       | K1       | O4  | $10 ~ \rm Mpc$     | for O4 high sensitivity simulations |
| `kagra_25Mpc.txt`       | K1       | O4  | $25 ~ \rm Mpc$     | for O4 top-end sensitivity simulations |
| `AplusDesign.txt`       | H1, L1   | O5  | $330 ~ \rm Mpc$    | A+ design target for O5 |
| `avirgo_O5high_NEW.txt` | V1       | O5  | $150 ~ \rm Mpc$    | target O5 sensitivity (high noise, low range) |
| `avirgo_O5low_NEW.txt`  | V1       | O5  | $260 ~ \rm Mpc$    | target O5 sensitivity (low noise, high range) |
| `kagra_80Mpc.txt`       | K1       | O5  | $80 ~ \rm Mpc$     | for O5 simulations |
| `kagra_128Mpc.txt`      | K1       | O5  | $128 ~ \rm Mpc$    | for O5 optimistic simulations |

### How to use

The pipeline only accepts file names with the following structure: `<detector-id>_<observing-run>.txt`. For any choice of file listed in the above table, please create a copy of the corresponding file in the folder where you plan to store the PSD files, and rename it to the corresponding structure.