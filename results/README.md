
# Results

## Results from testing different Sgr contamination removal strategies
### The Gaussian ball strategy is sufficient and not obviously biasing results at 40+ kpc distances.

All are run by following the pattern (after creating the `chains` directory and appropriate `prefixdir`)
```
python Reflex_fit_data.py processed_real/sgrtests/data_file.txt chains/prefixdir/
```

The chains can be inspected using `read_posterior.py`, which will produce the uncertainties.

All DR3 data (no Gaia designation) in the `processed_real/sgrtests` directory is created from the parent YPP24 data (the files without `dr2` in the string are DR3.). This data already has Sgr removed according to the PP21 Gaussian ball strategy. All DR2 data is created from the parent PP21 data, which does not have Sgr removed, and this we can test the effect of all Sgr cuts fairly. The isotropic halo is built from a standard EXP initial condition; there is nothing special about the halo itself. Any isotropic sphere of stars that sort of looks like the MW potential will do. This one isn't even that great (see note in the isotropic halo section).

### K giants, DR3 data

| Cut | Radius Limit | l | b | vtravel | vr | vphi | vtheta | sigmavlos | sigmul | sigmub | Nstars |
|-----|--------------|---|---|---------|----|------|--------|-----------|--------|--------|--------|
| 50 + data |
| YPP24all (quoted) | 50+ | 38 ± 11 | -37 ± 11 | 40 ± 7 | -9 ± 7 | -24 ± 6 | 17 ± 7 | 87 ± 4 | 70 ± 5 | 74 ± 5 | 340 |
| YPPKgiants (redo) | 50+ | 40 ± 11 | -36 ± 13 | 40 ± 7 | -5 ± 7 | -30 ± 6 | 19 ± 8 | 88 ± 5 | 61 ± 5 | 70 ± 5 | 253 | 
| B24| 50+ | 36 ± 23 | -45 ± 21 | 27 ± 9 | -1 ± 9 | -19 ± 7 | 17 ± 9 | 90 ± 5 | 58 ± 5 | 58 ± 5 | 160 |
| J21 | 50+ | 126 ± 26 | -18 ± 32 | 16 ± 6 | 2 ± 8 | 4 ± 6 | -16 ± 7 | 88 ± 5 | 36 ± 5 | 42 ± 5 | 153 |
| 40+ data | 
| YPP24all (quoted) | 40-50 | 
| YPPKgiants | 40+ | 61 ± 15 | -50 ± 12 | 30 ± 6 | -15 ± 5 | -23 ± 5 | 26 ± 6 | 94 ± 3 | 77 ± 4 | 84 ± 4 | 531 |
| B24 | 40+ | 78 ± 22 | -23 ± 25 | 19 ± 7 | -3 ± 7 | -10 ± 7 | 10 ± 7 | 97 ± 4 | 74 ± 4 | 62 ± 4 | 319 |
| J21 | 40+ | 141 ± 27 | 22 ± 23 | 13 ± 6 | -0.2 ± 6 | 10 ± 4 | -22 ± 5 | 92 ± 4 | 50 ± 3 | 56 ± 3 | 354 |

Summary:
1. All cuts return statistically significant detections of reflex motion.
2. The original PP21 Gaussian sphere cut and the spatial cut from B24 give consistent results. The J21 angular momentum cut is significantly different (and not in the direction of the C24 results?).

*This is appropos of nothing, but I was curious and had some time to kill while these runs executed. When running pyMultiNest, fully charged M1 mac power draw jumps from 7W to 22W (+15W). Models take approximately 10 minutes to run, or 2.5Wh=0.0025kWh. With the UK mix of power at 23.5p/kWh, the approximate cost is then  0.05p (that is, 20 calculations cost 1p). Assuming 100g/kWh, each calculation is about 0.2g CO2. According to Activity Monitor, I am successfully using all four performance cores to complete this calculation.*

### K giants, DR2 data


| Cut | Radius Limit | l | b | vtravel | vr | vphi | vtheta | sigmavlos | sigmul | sigmub | Nstars |
|-----|--------------|---|---|---------|----|------|--------|-----------|--------|--------|--------|
| 50+ data |
| PP21 | 50+ | 42 ± 10 | -30 ± 10 | 47 ± 7 | -8 ± 7 | -31 ± 8 | 10 ± 8 | 88 ± 4 | 66 ± 7 | 55 ± 6 | 263 
| none | 50+ | 34 ± 8 | -30 ± 10 | 49 ± 7 | -11 ± 7 | -29 ± 7 | 10 ± 8 | 87 ± 4 | 60 ± 7 | 53 ± 6 | 284
| B24 | 50+ | 38 ± 13 | -28 ± 14 | 41 ± 9 | -7 ± 8 | -24 ± 9 | 8 ± 9 | 90 ± 6 | 60 ± 9 | 48 ± 6 | 171
| J21 | 50+ | 115 ± 22 | 16 ± 21 | 21 ± 9 | 19 ± 9 | 10 ± 7 | -21 ± 8 | 86 ± 5 | 36 ± 8 | 32 ± 6 | 165
| 40+ data |
| PP21 reported| 40+ | 53 ± 9 | -28 ± 10 | 35 ± 5 | -17 ± 5 | -22 ± 5 | 18 ± 7 | 94 ± 3 | 81 ± 4 | 71 ± 3 | 543
| PP21 redo | 40+ | 51 ± 10 | -34 ± 11 | 35 ± 6 | -16 ± 5 | -26 ± 6 | 15 ± 6 | 94 ± 3 | 80 ± 5 | 73 ± 4 | 558
| none | 40+ | 31 ± 9 | -40 ± 8 | 41 ± 5 | -26 ± 5 | -29 ± 4 | 21 ± 6 | 92 ± 3 | 69 ± 4 | 77 ± 4 | 677
| B24 | 40+ | 57 ± 12 | -4 ± 14 | 32 ± 8 | -9 ± 7 | -13 ± 8 | -2 ± 8 | 98 ± 4 | 72 ± 6 | 54 ± 5 | 343
| J21 | 40+ | 129 ± 23 | 31 ± 18 | 17 ± 7 | 7 ± 6 | 14 ± 6 | -26 ± 6 | 92 ± 4 | 52 ± 4 | 49 ± 4 | 359

Summary:
1. DR2 vs DR3 is not a significant difference in any of the Sgr cut scenarios (PP21 and B24 agree; J21 does not. PP21 is even consistent-ish with not making an Sgr cut!).
2. In general, there isn't _that_ much Sgr-consistent material beyond 40 kpc, and particularly not past 50 kpc. The B24 cut could be combined with a very un-restrictive angular momentum cut (a larger PP21 ball enclosing the entirety of the Vasiliev 2021 model) and give a similar number to the PP21 cut at both 40+ and 50+. The J21 cut does remove a lot of stars, but I'm exploring that next...
3. I attribute the small change (statistically consistent) between the reported PP21 and the redone PP21 results to come from a smaller number of MultiNest walkers used in this analysis and slightly updated peculiar velocities.

### Mock isotropic halo test

From an isotropic MW-like halo, pull the first X number of halo particles in thick shells measured from the galactic centre. Then, subselect from that sample according to proposed Sgr cuts (currently just J21) to test the effect of different spatial and velocity cuts.

The J21 cut (angular momentum cut) eliminates approximately 1/3rd of the halo. It seems plausible this could create some bias, or at least change the results.

| Cut | Radius Limit | l | b | vtravel | vr | vphi | vtheta | sigmavlos | sigmul | sigmub | Nstars |
|-----|--------------|---|---|---------|----|------|--------|-----------|--------|--------|--------|
| 50+ data |
| none | 50-55 | 6 ± 101 | -1 ± 47 | 8 ± 9 | 2 ± 7 | 1 ± 16 | 2 ± 16 | 119 ± 6 | 15 ± 24 | 15 ± 14 | 253
| J21  | 50-55 | -14 ± 117 | -10 ± 50 | 10 ± 11 | 5 ± 9 | 10 ± 20 | -15 ± 21 | 120 ± 7 | 14 ± 16 | 14 ± 11 | 168
| J21 | 50-55 | -89 ± 69 | -23 ± 43 | 11 ± 9 | 0 ± 6 | 6 ± 13 | -31 ± 15 | 116 ± 5 | 26 ± 80 | 40 ± 73 | 315
| J21 | 50-55 | -94 ± 72 | -3 ± 41 | 7 ± 7 | 2 ± 5 | 6 ± 10 | -27 ± 11 | 118 ± 3 | 80 ± 32 | 61 ± 45 | 623
| 65+ data |
| none | 65-70 | -47 ± 80 | -18 ± 44 | 6 ± 6 | 9 ± 4 | 8 ± 11 | 17 ± 12 | 117 ± 3 | 100 ± 36 | 19 ± 59 | 811
| J21 | 65-70 | -53 ± 59 | -49 ± 21 | 16 ± 9 | 3 ± 5 | 4 ± 13 | 1 ± 15 | 116 ± 4 | 15 ± 15 | 17 ± 43 | 529
| 75+ data |
| none | 75-80 | -75 ± 36 | -16 ± 31 | 13 ± 8 | 4 ± 4 | 17 ± 13 | -5 ± 13 | 118 ± 3 | 16 ± 30 | 14 ± 14 | 807
| J21 | 75-80 | -73 ± 20 | -33 ± 19 | 29 ± 9 | 9 ± 6 | 51 ± 18 | 0 ± 19 | 114 ± 4 | 15 ± 14 | 14 ± 13 | 511
| 80+ data |
| none | 80-85 | -128 ± 43 | -47 ± 25 | 17 ± 8 | -2 ± 4 | 8 ± 15 | 13 ± 16 | 117 ± 3 | 14 ± 13 | 14 ± 11 | 787 
| J21 | 80-85 | -95 ± 55 | -60 ± 21 | 24 ± 10 | -3 ± 5 | 25 ± 18 | 23 ± 20 | 116 ± 4 | 14 ± 10 | 14 ± 11 | 534

Summary:
1. In models with no reflex motion, we largely recover no reflex motion!
2. Some evidence for a preferred direction when considering larger sample numbers, with travel signal increasing with distance. This seems bad! At the very least this test seems to demonstrate that a broad angular momentum cut changes the reflex detection in an aphysical way.
3. I note that my mocks are not extremely realistic; I used a single value for the uncertainties in velocities, and it looks as though I overestimated the proper motion uncertainties pretty badly (small hyperparameters = lots of width in the input distributions).

### Conclusions

It isn't obvious that the PP21 Sgr selection has changed any results. It does seem more true than ever now that the apex locatoin is difficult to pin down, and different selections may in fact make a difference at a larger-than-expected level. Future effort should be put in both continuing to stress test the method, while also characterising substructure via other means (chemical tagging?).