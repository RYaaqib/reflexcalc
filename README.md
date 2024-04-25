# reflexcalc

This reposity hosts the results of _Yaaqib,Petersen and Penarrubia (2024)_(https://arxiv.org/abs/2402.10841). We provide the code to reporduce the figures from the results (excluding the simulation lines, please get in touch with either me (Rashid Yaaqib) or Mike Petersen if you would like to use the reflex motion values calculated from simulations). 

This repo provides the tools to calculate a 'correction' arising from reflex motion for any given halo star with known 6D information. A file list is outlined below. 

We also provide the code to fit your own data using our same model, the fitting code is provided in full, but code to analyse the posterior chains will be added soon (under progress).

If you use this repo/codes for your reserach, please consider citing our publication at https://arxiv.org/abs/2402.10841 or the original paper _Petersen & Penarrubia_(2021, https://doi.org/10.1038/s41550-020-01254-3) (to be updated).

# File list
- coord.py
    - A file containing a list of coordinate conversions.
- Figures_paper.ipynb
    - A jupyter notebook that reporduces figures in the paper.
- calc_my_reflex.ipynb
    - A jupyter notebook that provides the steps to calculate the reflex motions for a given
    set of input data.
- genreflex.py 
    - A python file containing modified versions of the reflex motion model in reflex_fit_data.py, which 
    was modified to facilitate the plotting of the on-sky velocity maps
- Reflex_fit_data.py
    - A file to fit the reflex motion model given some input data. Usage:

         -In a terminal run 
         
         <code> python Reflex_fit_data.py /path/to/input.txt prefix</code>

        <code> /path/to/input.txt/ </code>: This is the relative path to the dat. The input data file should have the format described in the <code> Reflex_fit_data.py </code> docstring.

        <code> prefix </code>: This is the prefix for the __pyMultinest__ output files. Recomennded to be saved in a chains folder with chains/prefix.










