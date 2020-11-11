# cancer_covid_breast

Projections of impact of covid disruptions on breast cancer deaths in Scotland.
Accounts for differences in screen-detected and symptomatic breast cancers.

This code was adapted from code written to perform the analyses presented in the manuscript:

An inverse stage-shift model to estimate the excess mortality and health economic impact of delayed access to cancer services due to the COVID-19 pandemic

by Koen Degeling, Nancy N Baxter, Jon Emery, Fanny Franchini, Peter Gibbs, G Bruce Mann, Grant McArthur, Benjamin J Solomon, Maarten J IJzerman

For more information or when using/adapting the code, please refer to/cite the pre-print version of the manuscript available from MedRxiv: https://www.medrxiv.org/content/10.1101/2020.05.30.20117630v1

This code was adapted by Ewan Gray for the purpose of making predictions of cancer deaths in Scotland resuling from Covid-19 disruption to serivces leading to delays in diagnosis & treatment.

The model adaptations include:

    Scotland incidence and survival data as inputs
    TTI HR can be applied additionally to staging changes as an option.
    Seperation and recombination of results from sub-groups such as screen-detected and symptomatic in breast cancer.
    Selection of exponential rather than Gompertz as default option for projection of survival.
    Stage-shift from I to III/IV and II to III/IV are included. Shift to III and IV occur in proportion to the initial distribution of stage III/IV incident cases.
    Model is implemented as a function to allow convenient sensitivity analysis.
    A probabilistic sensivity analysis has been added.
