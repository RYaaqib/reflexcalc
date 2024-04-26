# reflexcalc

This reposity hosts the results of _Yaaqib,Petersen and Penarrubia (2024)_(https://arxiv.org/abs/2402.10841). We provide the code to reporduce the figures from the results (excluding the simulation lines, please get in touch with either me (Rashid Yaaqib) or Mike Petersen if you would like to use the reflex motion values calculated from simulations). 

This repo provides the tools to calculate a 'correction' arising from reflex motion for any given halo star with known 6D information. A file list is outlined below. 

We also provide the code to fit your own data using our same model, the fitting code is provided in full.

If you use this repo/codes for your reserach, please consider citing our publication at https://arxiv.org/abs/2402.10841 or the original paper _Petersen & Penarrubia_(2021, https://doi.org/10.1038/s41550-020-01254-3) (to be updated).

If there are any questions, I can be reached at <url>rashid.yaaqib@ed.ac.uk</url>.
---
# Acknowledgements

This repo uses the packages/repositories:
- [_pyMultinest_](https://github.com/JohannesBuchner/PyMultiNest): Buchner et al. [2014](http://www.aanda.org/articles/aa/abs/2014/04/aa22971-13/aa22971-13.html)
- [_ReflexMotion_](https://github.com/michael-petersen/ReflexMotion) Petersen & Penarrubia ([2021](https://ui.adsabs.harvard.edu/abs/2021NatAs...5..251P/abstract))
- _Numpy_ : Harris et. al ([2020](https://www.nature.com/articles/s41586-020-2649-2))
-_Matplotlib_: Hunter ([2007](https://ieeexplore.ieee.org/document/4160265))
---
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

        <code> /path/to/input.txt </code>: This is the relative path to the data. The input data file should have the format described in the <code> Reflex_fit_data.py </code> docstring.

        <code> prefix </code>: This is the prefix for the __pyMultinest__ output files. Recomennded to be saved in a chains folder with chains/prefix.
- read_posterior.py: A file containing helper functions to read the posterior chains returned by multinest. Technical note: Due to the wrapping of the $(\ell,b)_{\rm apex}$ parameters, the computation of the percentiles (width of the posteriors) needs a bit more work than using <code> np.percentile </code>. The percentiles for these quantities are done by shifting the posterior by the median such that it is centred at ~ 0. We only get the widths for the posterior from this, not the median.

---
# Required python packages.

As this is not a "package" but just some code on github - You will need the following packages:

- Numpy
- Scipy
- Pymultinest(only if you are using <code> Reflex_fit_data.py </code>).








