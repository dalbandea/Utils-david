module Utils

# import BDIO, FastGaussQuadrature
import ADerrors, Printf, Plots, LeastSquaresOptim, LaTeXStrings

include("UtilsFunctions.jl")
export UtilsFunc, MixedUtilsFunc, nparameters, FuncSeries, FuncPade, Correlator, SymCorrelator, ConstantFunc, Einstein, MixedCorrelator, ProjSymCorrelator, ProjMixedCorrelator, AutCorrScaling, AutCorrScalingwoOffset

include("UtilsFit.jl")
export fit_data, fit_data_yerr, fit_defs_yerr

include("UtilsPlot.jl")
export plot_fit, plot_data, plot_func, plot_fit!, plot_func!

#include("UtilsMSbar.jl")
#export betaMS

#include("UtilsInt.jl")
#export int_quad

##
## NOTE: A really useful module should allow for
##       different conventions beta(g), beta(alpha), ...
##
#include("UtilsBeta.jl")
#export GetSigma

##include("UtilsFlow.jl")
##export 



end # module
