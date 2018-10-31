Replication files for "The Power of Forward Guidance Revisited"
by Alisdair McKay, Emi Nakamura, and Jon Steinsson

Friday, January 29, 2016  (Final submission to AER)


Step 1: do the calculations by running main.m 
        (This takes about an hour because of figures 5 and 6.  If you want something faster see below.)
Step 2: Make the figures
        Figures 3 and 4 show the response to a 20-quarter-ahead 50bps cut.  To produce these figures run Figures_3_and_4.m
        Figures 5 and 6 show the response to a 50bps cut at different horizons. To produce these figures run Figures_5_and_6.m
        Figures 7, 8, and 9 show the ZLB episode.  To produce these figures: Figures_7_8_and_9.m

Step 3: Table 2 shows results for alternative calibrations.  
        edit FullModel/parameters/initparams.m (set options at top to select cases)
        run FullModel(false) to evaluate that case.



main.m solves for the transition path in response to interest rate changes at many horizons
in order to construct figures 5 and 6.  If you want to speed things up and omit these calculations
edit main.m and change FullModel(true) to FullModel(false).