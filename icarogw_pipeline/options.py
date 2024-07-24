
usage = """

\nWelcome to the ICAROGW runner helper!\n

    # ----- #
    # input #
    # ----- #

        output                      Directory from which input samples are read. The path is relative to the 'samples' namib directory. Default: ''
        
        injections-path             Flag to activate a hard comparison in case only two options are compared. This parameter is used only in violin and ridgeline plots. Default: 0
        injections-number           Flag to read the evidence from the samples. Default: 0
        snr-cut                     Flag to plot the prior. Option currently implemented only in corner plot and not fully implemented. Default: 0
        ifar-cut                    List of true values of the selected parameters. Currently implemented only with corner-sns=0. Default: []
        selection-effects-cut       Flag to read the evidence from the samples. Default: 0

        O3-cosmology                List of modes for which the QNMs are computed. This option is used only when QNMs are computed from {Mf, af}. Default: [(2,2)]
        simulation                  Flag to convert the damping time in [ms] and scale amplitudes as [1e-21]. The option is used to compare samples from Damped Sinusoids with other models. Default: 0
        data-path                   Flag to use pyRing fits to compute the QNMs. Default: 1

    # ----- #
    # model #
    # ----- #

        model-primary               Options: {'PowerLaw', 'PowerLaw-Gaussian', 'PowerLaw-GaussianRedshiftLinear', 'PowerLawRedshiftLinear-GaussianRedshiftLinear'}
        model-secondary             Options: {'PowerLaw', 'PowerLaw-Gaussian', 'MassRatio'}
        model-rate                  Options: {'MadauDickinson', 'BetaDistribution', 'BetaDistribution-Line', 'MadauDickinson-GammaDistribution', 'PowerLaw'}
        redshift-transition         Options: {'sigmoid', 'double-sigmoid', 'linear', 'linear-sinusoid'}

        positive-peak               List of parameters bounds that are used in the plots. Default: []. Syntax: [[a, b], [c, d], ...]
        low-smoothing               List of the selected parameters ordering. Default: []
        priors                      List of the compared parameters ordering. Default: []
        scale-free                  List of the compared parameters ordering. Default: []

    # ------- #
    # sampler #
    # ------- #

        nparallel                   Flag to save the imput samples filtered on the selected parameters. They are saved in 'output/reduced_posteriors'. Default: 0
        neffPE                      Flag to save the medians of the selected parameters. They are saved in 'output/output_medians'. Default: 0
        neffINJ                     Option to downsample the input samples, taking the corresponding percent value (downsampling=1, 100% of initial samples). Default: 1

        sampler                     Flag to read the evidence from the samples. Default: 0
        nlive                       Flag to read the evidence from the samples. Default: 0
        npool                       Flag to read the evidence from the samples. Default: 0

"""