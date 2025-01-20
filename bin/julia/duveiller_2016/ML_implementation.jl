### Dependencies

using Pkg

Pkg.activate("/Net/Groups/BGI/work_3/OEMC/oemc_sif/bin")

using LinearAlgebra, Optim, Plots, Dates

using DimensionalData, YAXArrays, Zarr, Statistics, LineSearches, Revise

using Rasters: Center
using Rasters

using MLJ, DataFrames, Tidier, StatsPlots, PlotlyJS

using ProgressMeter
# data at low resolution

lst_cube_low = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/lst_cube_low_july.zarr"))

ndwi_cube_low = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/ndwi_cube_low_july.zarr"))


ogvi_cube_low = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/ogvi_cube_low_july.zarr"))


sif_cube_low_july = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/sif_cube_low_july.zarr"))



# creating a dataframe with the information

low_res_df = CubeTable(lst = lst_cube_low, ndwi = ndwi_cube_low, ogvi = ogvi_cube_low, sif = sif_cube_low_july)

final_df = DataFrame(low_res_df[1])

# converting to longer format to filter easier

df_post = @chain final_df begin
    @pivot_longer(lst:sif, names_to = "variable", values_to = "value")
    @filter(!isnan(value))
    @pivot_wider(names_from = variable, values_from = value)
    @drop_missing()
end


@df df_post corrplot(cols(3:6))

features = [:lst, :ndwi, :ogvi, :sif]

# pair plot

PlotlyJS.plot(df_post, dimensions=features, kind="splom", marker=attr(opacity = 0.3, color = "black"))


sif, predictors = unpack(df_post[:,3:6], ==(:sif); rng=123);

# Checking which ML models we can use

models(matching(predictors,sif))

# I will focus XGBoostRegressor as a first approach

doc("XGBoostRegressor", pkg="XGBoost")
info("XGBoostRegressor", pkg="XGBoost")

XGBoostRegressor = @load XGBoostRegressor pkg=XGBoost

xgb = XGBoostRegressor()

evaluate(xgb, predictors, sif,
resampling=CV(nfolds = 10),
measure=[RootMeanSquaredError(), LPLoss()])

# Tuning model hyperparameters

xgb_m = EnsembleModel(model = xgb)

#=
from kaggle (https://www.kaggle.com/code/prashant111/a-guide-on-xgboost-hyperparameters-tuning/notebook)

 space={'max_depth': hp.quniform("max_depth", 3, 18, 1),
        'gamma': hp.uniform ('gamma', 1,9),
        'reg_alpha' : hp.quniform('reg_alpha', 40,180,1),
        'reg_lambda' : hp.uniform('reg_lambda', 0,1),
        'colsample_bytree' : hp.uniform('colsample_bytree', 0.5,1),
        'min_child_weight' : hp.quniform('min_child_weight', 0, 10, 1),
        'n_estimators': 180,
        'seed': 0
    }
=#

# grid

r1 = range(xgb_m, :(model.max_depth), lower=3, upper=18);
r2 = range(xgb_m, :(model.gamma), lower=1, upper=9);
r3 = range(xgb_m, :(model.alpha), lower=40, upper=180);
r4 = range(xgb_m, :(model.lambda), lower=0, upper=1);
r5 = range(xgb_m, :(model.colsample_bytree), lower = 0.5, upper =1);
r6 = range(xgb_m, :(model.min_child_weight), lower = 0, upper =10);


self_tuning_xgbost = TunedModel(
    model=xgb_m,
    tuning=Grid(goal=30),
    resampling=CV(nfolds=10),
    range=[r1, r2, r3, r4, r5, r6],
    measure=RootMeanSquaredError());

mach = machine(self_tuning_xgbost, predictors, sif);
#fit!(mach, verbosity=1);

Plots.plot(mach)

report(mach)

#MLJ.save("/Net/Groups/BGI/work_3/OEMC/oemc_sif/results/model_xgboost_v1.jls", mach.fitresult)


mach = machine("/Net/Groups/BGI/work_3/OEMC/oemc_sif/results/model_xgboost_v1.jls")

####### Downscaling SIF #############


lst_cube_high_july = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/lst_cube_high_july.zarr"))

ogvi_cube_high_july = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/ogvi_cube_high_july.zarr"))

ndwi_cube_high_july = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/ndwi_cube_high_july.zarr"))


predict(mach, DataFrame(lst=200, ndwi=0.3, ogvi=0.5))

#=
indims = (InDims(),InDims(),InDims())
outdims = (OutDims())

function sif_downscaling(xout, lst, ndwi, ogvi; model)

    if all(!ismissing, [lst[1], ogvi[1], ndwi[1]]) && all(!isnan, [lst[1], ogvi[1], ndwi[1]])
        xout .= predict(model, DataFrame(lst=lst[1], ndwi = ndwi[1], ogvi = ogvi[1]))
    end
end

mapCube(sif_downscaling, (lst_cube_high_july, ndwi_cube_high_july, ogvi_cube_high_july), indims = indims, outdims = outdims; model = mach)
=#


df_high_res = CubeTable(lst = lst_cube_high_july, ndwi = ndwi_cube_high_july, ogvi = ogvi_cube_high_july)

df_high_res_f = DataFrame(df_high_res[1])

sif_pred = predict(mach, df_high_res_f[:,1:3])

test = reshape(sif_pred, (4800,4800))

test2 = replace(test, test[1,1] => NaN)

p1 = Plots.heatmap(test2, clim = (0,2.8), title = "SIF (Sentinel-5p) Downscaled 1km (2018-07)")

p2 = Plots.heatmap(sif_cube_low_july.data[:,:], clim = (0,2.8), title = "SIF (Sentinel-5p) 10 km (2018-07)")

Plots.plot(p2,p1)

Plots.plot!(size=(1000,400))

Plots.savefig("/Net/Groups/BGI/work_3/OEMC/oemc_sif/results/sif_downscaling_XGBoost.png")



#######################################################


# MLJ tutorial

iris = load_iris();

selectrows(iris, 1:3)  |> pretty

import DataFrames
iris = DataFrames.DataFrame(iris);

ytmp, X_tmp = unpack(iris, ==(:target); rng=123);

X_tmp

models(matching(X_tmp,ytmp))

doc("DecisionTreeClassifier", pkg="DecisionTree")
info("DecisionTreeClassifier", pkg="DecisionTree")
Tree = @load DecisionTreeClassifier pkg=DecisionTree

tree = Tree()


evaluate(tree, X_tmp, ytmp,
resampling=CV(shuffle=true),
measures=[log_loss, accuracy],
verbosity=0)


typeof(y)

target_scitype(tree)

scitype(y)


yint = int.(ytmp)


scitype(yint)


mach = machine(tree, X_tmp, ytmp)

train, test = partition(eachindex(ytmp), 0.7); # 70:30 split

MLJ.fit!(mach, rows=train)

yhat = MLJ.predict(mach, X_tmp[test,:])

yhat[3:5]

log_loss(yhat, y[test])

# predicted probabilities of virginica

broadcast(pdf, yhat[3:5], "virginica")

# predicted probability of observed class

broadcast(pdf, yhat, y[test])[3:5]


mode.(yhat[3:5])


predict_mode(mach, X_tmp[test[3:5],:])

evaluate!(mach, resampling=Holdout(fraction_train=0.7),
measures=[log_loss, accuracy],
verbosity=0)

tree.max_depth = 3

evaluate!(mach, resampling=Holdout(fraction_train=0.7),
measures=[log_loss, accuracy],
verbosity=0)

###############################################


# crea