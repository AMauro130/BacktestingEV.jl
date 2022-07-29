module BacktestingEV

export back_real_vs_meanPred_plot, back_real_vs_meanpred
export back_periodicity
export back_fitting_ev, back_perc_ev
export back_testingEV

include("back_testing_fct.jl")
include("utils.jl")
include("criterias.jl")
include("threshold_selection.jl")

end
