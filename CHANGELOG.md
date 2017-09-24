CHANGELOG: fastspt, a library fitting realistic jump length distribution to SPT data
------------------------------------------------------------------------------------

## v12.0: Fix some parameters in the model -- (2017-09-24)
This version allows to fix diffusion coefficients in the model (in case
those parameters have already been estimated from an independent source.

To use it, simply set the lower bound and the upper bound of the corresponding 
diffusion coefficient to the independently estimated value.

**Implementation details:** Look at the difference between LB and UB, if 
under a float threshold, consider that this parameter should be fixed, and 
declare it as is.

**API change:** this should allow not to break the API, and provide natural
evolution for aplications that depend on fastspt, such as Spot-On.



## v7.?: Plotting functions (2017-03-19)
- Started a Python port of Anders' code
- Including the plotting system


## v7: UPDATE Python port (2017-03-12)
- Started a port in Python


## v6: UPDATE June 13, 2016 (2016-06-13)
Performed control experiments to measure the axial detection slice that
can actually be determined using fixed cells labeled with JF549 or
JF646 and z-stack imaging with 20 nm z-stack. It is difficult to say
with absolue certainty, but dZ appears to be around 700 nm. Performed
new Monte Carlo simulations and determined a new empirical
relationship assuming dT = 4.477 ms and 1 gap allowed:

```{matlab}
dZ_corr = dZ + a*sqrt(D) + b
dZ_corr = 700 nm + 0.15716*sqrt(D) + 208.11 nm
```

Updated the code to reflect this. `SS_2State_model_Z_corr_v4` reflect
this update. 


## v5: UPDATE June 10, 2016 (2016-06-10)
Re-wrote data loading section to allow for loading in data from
multiple different experiments and fitting a curve to all of the data.
This is mainly for figure visualization purposes, so do not use
Single-Cell fitting when doing this. 

## v4: UPDATE June 3, 2016 (2016-06-03)
Performed Monte Carlo simulations to find the optimal Z-corr for your
conditions: dZ = 0.8 um; dT = 4.5 ms; 1 gap allowed;
So now using the new empirical relationships:
`dZ_corr = dZ + 0.144*sqrt(D) + 0.271;`


## v3: UPDATE June 2, 2016 (2016-06-02)
Modify the code to allow cell-by-cell analysis. This will help
identifying any biasing aspects of any particular movies and movies
where the particle density was too high or too low can be excluded. 


## v2: UPDATE May 8, 2016 (2016-05-08)
Since you are not accounting for transitions between bound and unbound
within DeltaTau because simulations showed this was negligble for CTCF,
there is no need to use $k_ON$ and $k_OFF$ when you do the fitting. Thus,
in "v2" we now use Fraction_Bound and (1 - Fraction_Bound) to do the fitting,
instead of $k_ON$ and $k_OFF$. This way, there are only three free parameters
and you can just calculate t_search for $k_OFF$ (which you know) and
Fraction bound. 


## UPDATE May 4, 2016 - Z_corr & PDF/CDF fitting (2016-05-04)
Correct for molecules that are moving in the axial dimension. The moving
fraction is more likely to move out of focus. Use the procedure of Kues
and Kubitscheck assuming an absorbing boundary and then use the empirical
relationship from Davide Mazza:

```{matlab}
DeltaZ_USE = DeltaZ_TRUE + 0.21 * sqrt(D); % (v4 uses different expression)
```

Also, you can now do both PDF and CDF fitting. The variable `PDF_or_CDF`
determines whether you are fitting to a histogram or to a CDF function:

- `PDF_or_CDF = 1` ==> perform histogram fitting
- `PDF_or_CDF = 2` ==> perform CDF fitting

	This function now includes a correction for the fact that free
	molecules are more likely to be lost and move out of focus. 

