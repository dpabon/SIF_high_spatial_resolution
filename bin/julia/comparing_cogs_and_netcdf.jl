using Pkg

Pkg.activate("/Net/Groups/BGI/work_3/OEMC/oemc_sif/bin/julia")

using DimensionalData
using YAXArrays
using ArchGDAL
using NetCDF
using Dates
using Plots
using DataFrames
using CSV
using CairoMakie
using Makie
using GeoMakie

cogpath  = ("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/OEMC_app_data/COG_SIF_files/")

files = readdir(cogpath)

cubes = Cube.(readdir(cogpath, join = true));

cubes

files[1][26:32]
files[1:48]
# extracting time from files

cog_dates = DateTime.(SubString.(files[:], 27,34), "yyyymmdd") 

cog_cube = cat(cubes..., dims=Ti(cog_dates));


lookup(cog_cube, :Ti)

nc_path = "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/OEMC_app_data/Monitor_23_World-carbon_flux_monitor_netcdf/"

files_nc = readdir(nc_path)

test = open_dataset(nc_path.*files_nc[5])
Cube(test[["SIF"]])



nc_cubes = Cube.([open_dataset(nc_path.*files_nc[i])[["SIF"]] for i in eachindex(files_nc)])

nc_cube = cat(nc_cubes..., dims=Ti)


lon = lookup(nc_cube, :lon)
lat = lookup(nc_cube, :lat)


data = nc_cube[Ti = At(DateTime("2010-04-19"))].data[:,:]

data_f = replace(data, data[1,1] => missing)

CairoMakie.activate!()

fig, ax, plt = CairoMakie.heatmap(data_f; colormap = :seaborn_icefire_gradient,
    axis = (; aspect=DataAspect()),
    figure = (; size = (1200,600), fontsize=24))
fig

fig = Figure(;size=(1200,600))
ax = GeoAxis(fig[1,1], dest = "+proj=latlon")
CairoMakie.heatmap!(ax, lon, lat, data_f; colormap = :seaborn_icefire_gradient, shading=false)
#cl=lines!(ax, GeoMakie.coastlines(), color = :white, linewidth=0.85)
#CairoMakie.translate!(cl, 0, 0, 1000)
fig

#################### Comparing Netcdf with COG files

dates_org = lookup(nc_cube, :Ti)

start_date = dates_org .- Day(8)

end_date = dates_org .+ Day(7)

out_csv = DataFrame(start_date = start_date, end_date = end_date)

CSV.write("/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/OEMC_app_data/start_end_dates.csv", out_csv)



#heatmap(nc_cube[Ti = 1].data)
#heatmap(cog_cube[Ti = 1].data)


nc_cube[Ti = 1].data[1,1]

nc_cube[Ti = 1].data[1:200,1:200]
heatmap(cog_cube[Ti = 1][1000:1800, 1500:1800].data)

all(ismissing, cog_cube[Ti = 1])



DateTime("1970-01-01") +Day(14718)