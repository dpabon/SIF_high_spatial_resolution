

### Dependencies


using Pkg

Pkg.activate("/Net/Groups/BGI/work_3/OEMC/oemc_sif/bin")

using LinearAlgebra, Optim, Plots, Dates

using DimensionalData, YAXArrays, Zarr, Statistics, LineSearches, Revise

using Rasters: Center
using Rasters


### Functions from Duveiller and Cescatti 2016 (DOI:10.1016/j.rse.2016.04.027)



function focalWeightFilter(x)
    scale_fact = x
    mWm = scale_fact * 3
    xCent = repeat([-scale_fact, 0, scale_fact], inner=(x,mWm))'
    yCent = repeat([scale_fact, 0, -scale_fact], inner=(x, mWm))
    vPos = collect(1:scale_fact) .- (scale_fact / 2 + 0.5)
    xPos = repeat(vPos, outer=(x,mWm))'
    yPos = repeat(sort(vPos, rev=true), outer = (x,mWm))
    dista = -(sqrt.((xPos - xCent).^2 + (yPos - yCent).^2) ./ scale_fact)
    dista[isinf.(dista)] .= 0
    dista = dista .- (minimum(dista) - 0.000001)
    dista = dista .+ 0.1
    weights_matrix = fill(NaN, (mWm, mWm))
    
    for n in 1:scale_fact
        j = n - ((n - 1) ÷ scale_fact) * scale_fact
        k = scale_fact .* (0:2) .+ j
        
        for m in 1:scale_fact
            o = scale_fact .* (0:2) .+ m
            distaSum = sum(dista[o, k])
            weights_matrix[o, k] .= dista[o, k] ./ distaSum
        end
    end
    
    return weights_matrix
end

focalWeightFilter(3)

# Weighted focal filter (gaussian)

function gaussianWeightFilter(x)
    scale_fact = x
    mWm = scale_fact * 3
    yM = repeat(1:mWm, outer=(1,mWm))
    xM = repeat(1:mWm, outer=(1,mWm))'
    
    dista = exp.(-(((xM .- (scale_fact * 1.5)).^2) .+ ((yM .- (scale_fact * 1.5)).^2))./ (2 * scale_fact^2))

    weights_matrix = fill(NaN, mWm, mWm)

    for n in 1:scale_fact
        j = n - ((n - 1) ÷ scale_fact) * scale_fact
        k = scale_fact .* (0:2) .+ j
        for m in 1:scale_fact
            o = scale_fact .* (0:2) .+ m
            distaSum = sum(dista[o, k])
            weights_matrix[o, k] = dista[o, k] / distaSum
        end
    end

    return weights_matrix
end


# cost function

function costfun(b,VI,ET,LST,sif_ob)

    sif_pred = vegetation(VI, b[1], b[2]) .* water(ET, b[3], b[4]) .* temperature(LST, b[5], b[6])
    return sum((sif_pred .- sif_ob).^2)
    
end


# model components

function vegetation(VI, b1, b2)
    b2 .* (VI.^b1)
end

function water(ET, b3, b4)
    1 ./(1 .+exp.(b3.*(-ET.+b4)))
end

function temperature(LST, b5, b6)
    exp.(-0.5.*((LST.+b5)./b6).^2) 
end

function sif(VI, ET, LST, b)
    return vegetation(VI, b[1], b[2]) .* water(ET, b[3], b[4]) .* temperature(LST, b[5], b[6])
end

## Testing the functions ##

param_ini=[1, 2, 50.0, 0, -295, 10]

param_min=[0.5, 0.1, 0.0, -1, -310, 1]
param_max=[1.5, 5, 500.0, 1, -290, 50]


vi = [0.3,0.2,0.4, 0.8, 0.2]
et = [100.,200.,120., 40., 50.]
lst  = [28.,30.,25., 21., 28.]

sif_ob = sif(vi, et, lst, [1.35, 3, 120, -0.5, -291, 20])


inner_optimizer = LBFGS()

test = optimize(b -> costfun(b, vi, et, lst, sif_ob), param_min, param_max, param_ini, Fminbox(inner_optimizer))


Optim.minimizer(test)

Plots.scatter([sif_ob sif(vi,et,lst,Optim.minimizer(test))])



## Preprocessing

### Gridding TROPOMI data


#sif_gridding = `nohup julia ./gridL2_Dates.jl --latMin 29 --latMax 71 --lonMin -30 --lonMax 89 --dLat 0.05 --dLon 0.05  --dDays 8 --startDate 2020-01-01 --stopDate 2020-12-31  --Dict specifications.json  -o /Net/Groups/BGI/scratch/jgens/TROPOSIF/gridded/Europe_h17-h20_v2-v5/Europe_h17-h20_v2-v5_005_8daily_2020.nc &`

# run the gridding of TROPOMI product.

# run(sif_gridding)


# See preprocessing_jhoannes files for the migrations steps from Netcdf to Zarr


## Data
#=
In the following example we will downscale TROPOSIF using Sentinel-3 satellite data and MODIS for a month (2018-07).
=#
### SIF low resolution


sif_cube = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/SIF_TROPOMI_gridded_10km_filtered_europe_2018-05-04_2020-12-26.zarr/"))


### selecting only SIF_743 and 2018-07


sif_cube_low = sif_cube[Ti = Between(DateTime("2018-07-01"),DateTime("2018-07-31")), Variable = At("SIF_743")]

heatmap(sif_cube_low.data[:,:,1])


### LST high resolution


lst_cube = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/LST_Sentinel3_gridded_1km_europe_2018-01-04_2020-12-29.zarr/"))

## selecting only 2018-07

lst_cube_high = lst_cube[Ti = Between(DateTime("2018-07-01"),DateTime("2018-07-31")), Variable = At("LST")]


### (WATER) MODIS NDWI


modis_cube = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/MC43A4_gridded_1km_europe_2018-05-01_2020-12-31_YAXA.zarr/"))


## selecting only 2018-07


modis_cube_high = modis_cube[Ti = Between(DateTime("2018-07-01"),DateTime("2018-07-31"))]

# Computing NDWI

function ndwi(out, b4,b2)
    if !isnan.(b4) && !isnan.(b2) && !isinf.(b4) && !isinf.(b2)
        out .= (b4[:]-b2[:])/(b4[:] + b2[:])

    else
        out .= NaN32
    end
end

indims = (InDims(), InDims())

outdims = OutDims(Dim{:Variable}(["NDWI"]))

ndwi_cube_high = mapCube(ndwi, (modis_cube_high[Variable = At("4")], modis_cube_high[Variable = At("2")]), indims = indims, outdims=outdims)

# quick check
heatmap(ndwi_cube_high.data[:,:,3,1])


### (VEGETATION) Sentinel-3 OTCI (Terrestrial Chlorophyll Index)


sentinel_3_cube = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/OLCI_Sentinel3_gridded_1km_europe_2018-01-04_2020-12-29.zarr/"))

sentinel_3_cube.Variable

### selecting only OGVI and 2018-07

ogvi_cube_high = sentinel_3_cube[Variable = At("OGVI"), Ti = Between(DateTime("2018-07-01"),DateTime("2018-07-31"))]


#=
Collapsing time dimension. As data for each one of the product has a different length, and interpolation with 4 points is risky let's estimate the mean for the month of july
=#

lookup(sif_cube_low, :Ti)
lookup(lst_cube_high, :Ti)
lookup(ndwi_cube_high, :Ti)
lookup(otci_cube_high, :Ti)

function temporal_mean(out, in)
    if !all(isnan, in)
        out .= mean(filter(!isnan, in))
    else
        out .= NaN32 
    end
    
end

indims = InDims(:Ti)
outdims = OutDims()

sif_cube_low_july = mapCube(temporal_mean, sif_cube_low, indims = indims, outdims = outdims)



lst_cube_high_july = mapCube(temporal_mean, lst_cube_high, indims = indims, outdims = outdims)

ndwi_cube_high_july = mapCube(temporal_mean, ndwi_cube_high, indims = indims, outdims = outdims)

ogvi_cube_high_july = mapCube(temporal_mean, ogvi_cube_high, indims = indims, outdims = outdims)


### saving data data



savecube(sif_cube_low_july, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/sif_cube_low_july.zarr"; overwrite = true)
savecube(lst_cube_high_july, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/lst_cube_high_july.zarr"; overwrite = true)
savecube(ogvi_cube_high_july, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/ogvi_cube_high_july.zarr"; overwrite = true)
savecube(ndwi_cube_high_july, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/ndwi_cube_high_july.zarr"; overwrite = true)



### reading data


sif_cube_low_july = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/sif_cube_low_july.zarr"))

# let's match X and Y in the other cubes

#=
sif_cube_low_july = sif_cube_low_july[X=1:100, Y = 1:100]

axlist = (
    Dim{:X}(1:100),
    Dim{:Y}(1:100))

sif_cube_low_july = YAXArray(axlist, sif_cube_low_july.data[:,:])
=#

lst_cube_high_july = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/lst_cube_high_july.zarr"))

ogvi_cube_high_july = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/ogvi_cube_high_july.zarr"))

ndwi_cube_high_july = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/ndwi_cube_high_july.zarr"))


sif_small_area = sif_cube_low_july[X = 1:100, Y = 1:100]

lst_small_area = lst_cube_high_july[X = 1:100, Y = 1:100]

ogvi_small_area = ogvi_cube_high_july[X = 1:100, Y = 1:100]

ndwi_small_area = ndwi_cube_high_july[X = 1:100, Y = 1:100]

## Aggregating data to match 10 km SIF resolution

function aggregate_space(out, incube; scale_factor, axes)
    test = Raster(replace(incube, NaN32 => missing), axes)
    test_out = aggregate(mean, test, (X(scale_factor), Y(scale_factor)); skipmissingval=true)
    out .= replace(test_out.data, missing => NaN32)
end


indims = InDims(:X,:Y)

outdims = OutDims(Dim{:X}(1:480), Dim{:Y}(1:480))

axes = lst_cube_high.axes[1:2]

scale_factor = 10



### LST aggregation


lst_cube_low = mapCube(aggregate_space, lst_cube_high_july, indims = indims, outdims = outdims; scale_factor=scale_factor, axes=axes)
   
heatmap(lst_cube_low.data[:,:])

savecube(lst_cube_low, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/lst_cube_low_july.zarr")


lst_cube_low = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/lst_cube_low_july.zarr"))

### NDWI aggregation

ndwi_cube_low = mapCube(aggregate_space, ndwi_cube_high_july, indims = indims, outdims = outdims; scale_factor=scale_factor, axes=axes)

heatmap(ndwi_cube_low.data[:,:])

savecube(ndwi_cube_low, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/ndwi_cube_low_july.zarr")

ndwi_cube_low = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/ndwi_cube_low_july.zarr"))

### OGVI
ogvi_cube_low = mapCube(aggregate_space, ogvi_cube_high_july, indims = indims, outdims = outdims; scale_factor=scale_factor, axes=axes)

heatmap(ogvi_cube_low.data[:,:])

savecube(ogvi_cube_low, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/ogvi_cube_low_july.zarr")

ogvi_cube_low = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/low_res/ogvi_cube_low_july.zarr"))


function param_optim(out, sif, vegetation, water, temperature; param_ini = [1.,2.,50.,0.,-295.,10.], param_min = [0.5,0.1,0,-1.,-310.,1.], param_max = [1.5,5.,500.,1.,-290.,50.], min_obs = 21, inner_optimizer = LBFGS())
    
    
    # checking that we have at least min_obs data points for each variable in the local window
    
    if !all(isnan, sif) && !all(isnan, vegetation) && !all(isnan, water) && !all(isnan, temperature)
        if length(filter(!isnan, vec(sif))) >= min_obs && length(filter(!isnan, vec(vegetation))) >= min_obs && length(filter(!isnan, vec(water))) >= min_obs && length(filter(!isnan, vec(temperature))) >= min_obs
            
            if length(filter(!isnan, vec(sif))) == length(filter(!isnan, vec(vegetation))) && length(filter(!isnan, vec(sif)))== length(filter(!isnan, vec(water))) && length(filter(!isnan, vec(sif))) == length(filter(!isnan, vec(temperature)))
                
                vi = filter(!isnan, vec(vegetation))
                agua = filter(!isnan, vec(water))
                lst = filter(!isnan, vec(temperature))
                sif_ob = filter(!isnan, vec(sif))
                
                test = optimize(b -> costfun(b, vi, agua, lst, sif_ob), param_min, param_max, param_ini, Fminbox(inner_optimizer))
                
                out .= Optim.minimizer(test) 
                
                
            else
                out .= NaN32
                
            end
        else
            out .= NaN32
        end 
    else
        out .= NaN32
    end
    
end


window_edge = 5


# cambiar a mayor window size!!!

if isodd(window_edge) 
    pre_step = after_step = floor(window_edge / 2) 
else 
    pre_step = after_step = floor(window_edge / 2) - 1
end

indims = (InDims(MovingWindow(:X, pre_step, after_step), MovingWindow(:Y, pre_step, after_step), window_oob_value = NaN), InDims(MovingWindow(:X, pre_step, after_step), MovingWindow(:Y, pre_step, after_step), window_oob_value = NaN),InDims(MovingWindow(:X, pre_step, after_step), MovingWindow(:Y, pre_step, after_step), window_oob_value = NaN),InDims(MovingWindow(:X, pre_step, after_step), MovingWindow(:Y, pre_step, after_step), window_oob_value = NaN))

outdims = OutDims(Dim{:parameters_optim}(["b1", "b2", "b3", "b4", "b5", "b6"]))

parameters_cube = mapCube(param_optim, (sif_cube_low_july[X=1:100, Y = 1:100], ogvi_cube_low[X=1:100, Y = 1:100], ndwi_cube_low[X=1:100, Y = 1:100], lst_cube_low[X=1:100, Y = 1:100]), indims = indims, outdims = outdims, inner_optimizer = LBFGS()) 


heatmap(parameters_cube.data[1,:,:])

#= savecube(parameters_cube, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/parameters_cube_low_OGVI.zarr")

parameters_cube = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/parameters_cube_low.zarr")) =#

# disaggregate coefficients to match the high resolution

function disaggregate_space(out,incube; scale_factor, axes) 
    test = Raster(replace(incube, NaN32 => missing), axes) 
    test_out = disaggregate(Center(), test, (scale_factor, scale_factor)) 
    out .= replace(test_out.data, missing => NaN32)

end

indims = InDims(:X,:Y)

outdims = OutDims(Dim{:X}(1:4800), Dim{:Y}(1:4800))

axes = lst_cube_low.axes[1:2]

scale_factor = 10

# disaggregate parameters cube

parameters_cube_high = mapCube(disaggregate_space, parameters_cube, indims = indims, outdims = outdims; scale_factor=scale_factor, axes=axes)

heatmap(parameters_cube_high.data[:,:,1]) 

savecube(parameters_cube_high, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/parameters_cube_high_OGVI.zarr/"; overwrite = true) 
parameters_cube_high = Cube(open_dataset("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/parameters_cube_high.zarr/"))


# function to estimate sif using high resolution data

function sif_downscaling(out, vegetation, water, temperature, parameters) 
    if !all(isnan, vegetation) && !all(isnan, water) && !all(isnan, temperature) && !all(isnan, parameters)
        
        if length(filter(!isnan, vec(vegetation)))== length(filter(!isnan, vec(water))) && length(filter(!isnan, vec(vegetation))) == length(filter(!isnan, vec(temperature)))
            
            b_optim = [mean(parameters[:,:,1]),
            mean(parameters[:,:,2]),
            mean(parameters[:,:,3]),
            mean(parameters[:,:,4]),
            mean(parameters[:,:,5]),
            mean(parameters[:,:,6])]
            length(b_optim)
            
            out .= sif(mean(filter(!isnan, vegetation)), mean(filter(!isnan, water)), mean(filter(!isnan, temperature)), b_optim)
        else
            out .= NaN32
        end 
    else
        out .= NaN32
    end
end

# using another moving window of 3x3 to Downscale SIF

window_edge = 3

if isodd(window_edge) 
    pre_step = after_step = floor(window_edge / 2) 
else 
    pre_step = after_step = floor(window_edge / 2) - 1 
end

indims = (InDims(MovingWindow(:X, pre_step, after_step), MovingWindow(:Y, pre_step, after_step), window_oob_value = NaN), InDims(MovingWindow(:X, pre_step, after_step), MovingWindow(:Y, pre_step, after_step), window_oob_value = NaN),InDims(MovingWindow(:X, pre_step, after_step), MovingWindow(:Y, pre_step, after_step), window_oob_value = NaN),InDims(MovingWindow(:X, pre_step, after_step), MovingWindow(:Y, pre_step, after_step), window_oob_value = NaN, :parameters_optim))

outdims = OutDims(Dim{:SIF}(["sif_downscaled"]))

sif_cube_high = mapCube(sif_downscaling, (ogvi_cube_high_july, ndwi_cube_high_july, lst_cube_high_july, parameters_cube_high), indims = indims, outdims = outdims)

# low resolution

heatmap(sif_cube_low_july.data[:,:])

#max_value = maximum(filter(!isnan,sif_cube_low_july.data[:,:])) 

#min_value = minimum(filter(!isnan,sif_cube_low_july.data[:,:]))

# high resolution

heatmap(sif_cube_high.data[1,:,:])

# OTCI

heatmap(otci_cube_high_july.data[:,:])

# LST

heatmap(lst_cube_high_july.data[:,:])

# NDWI

heatmap(ndwi_cube_high_july.data[:,:])

#parameters

heatmap(parameters_cube.data[1,1:200,1:200]) 
heatmap(parameters_cube.data[2,:,:]) 
heatmap(parameters_cube.data[3,:,:]) 
heatmap(parameters_cube.data[4,:,:]) 
heatmap(parameters_cube.data[5,:,:]) 
heatmap(parameters_cube.data[6,:,:])