
# Results

## Results testing different Sgr contamination removal strategies

All are run by following the pattern (after creating the `chains` directory and appropriate `prefixdir`)
```
python Reflex_fit_data.py processed_real/sgrtests/data_file.txt chains/prefixdir/
```

### K giants, DR3 data

| Cut | Radius Limit | l | b | vtravel | vr | vphi | vtheta | sigmavlos | sigmul | sigmub | Nstars |
|-----|--------------|---|---|---------|----|------|--------|-----------|--------|--------|--------|
| YPP24all (quoted) | 50+ | 38 ± 11 | -37 ± 11 | 40 ± 7 | -9 ± 7 | -24 ± 6 | 17 ± 7 | 87 ± 4 | 70 ± 5 | 74 ± 5 | 340 |
| B24| 50+ | 36 ± 23 | -45 ± 21 | 27 ± 9 | -1 ± 9 | -19 ± 7 | 17 ± 9 | 90 ± 5 | 58 ± 5 | 58 ± 5 | 160 |
| J21 | 50+ | 126 ± 26 | -18 ± 32 | 16 ± 6 | 2 ± 8 | 4 ± 6 | -16 ± 7 | 88 ± 5 | 36 ± 5 | 42 ± 5 | 153 |
| no additional | 50+ | 40 ± 11 | -36 ± 13 | 40 ± 7 | -5 ± 7 | -30 ± 6 | 19 ± 8 | 88 ± 4 | 61 ± 5 | 69 ± 5 | 255 |
| YPPKgiants (redo) | 50+ | 253 | 
| B24 | 40+ | 78 ± 22 | -23 ± 25 | 19 ± 7 | -3 ± 7 | -10 ± 7 | 10 ± 7 | 97 ± 4 | 74 ± 4 | 62 ± 4 | 319 |
