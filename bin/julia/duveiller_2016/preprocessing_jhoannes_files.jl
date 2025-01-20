using Pkg

Pkg.activate("/Net/Groups/BGI/work_3/OEMC/oemc_sif/bin")

using Zarr, NetCDF, Glob, Dates, YAXArrays, DimensionalData, ProgressMeter


## Creating a zarr cube with the proper format and metadata for each one of the products 


## TROPOSIF gridded (10 km, 8 days)

sif_files = glob("SIF_*.nc", "/Net/Groups/BGI/data/DataStructureMDI/DATA/Incoming/TROPOSIF/gridded/Europe_h17-h20_v2-v5/SIF_10km_filtered/")

sif_files

# Extracting latitude and longitude vectors

grid = Cube(open_dataset("/Net/Groups/BGI/data/DataStructureMDI/DATA/Incoming/TROPOSIF/gridded/Europe_h17-h20_v2-v5/SIF_10km_filtered/grid_h17-h20_v2-v5_10km.nc"))

grid.data[:,:,1]


lat = grid.data[:,:,1]'

lon = grid.data[:,:,2]'

# time dimension

time = sif_files

time = DateTime.([time[i][end-21:end-12] for i in eachindex(time)])

# variable dimension

first_sif = Cube(open_dataset(sif_files[1]))

Variable = collect(lookup(first_sif, :Variable))

push!(Variable, "lon", "lat")

# YAXArrays dimension

data_empty = Array{Float32, 4}(undef, 480,480,length(time),5)

axlist = (
    Dim{:X}(1:480),
    Dim{:Y}(1:480),
    Dim{:time}(time),
    Dim{:Variable}(Variable))



for i in eachindex(sif_files)
    tmp_cube = Cube(open_dataset(sif_files[i]))
    new_cube = collect(tmp_cube.data)

    for d in 1:3
        new_cube[:,:,d] = reverse(new_cube[:,:,d], dims=2)'
    end

    data_empty[:,:,i,1:3] = new_cube

    data_empty[:,:,i,4] = lon

    data_empty[:,:,i,5] = lat

end


props = Dict(
    "X" => "lon",
    "Y" => "lat",
    "time" => "days",
    "lon" => "longitude",
    "lat" => "latitude",
    "gridding" => "10km",
    "proj" => "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
    "Author" => "MPI-BGC, Daniel Pabon (dpabon@bgc-jena.mpg.de)"
)

ds = YAXArray(axlist, data_empty, props)

ds = setchunks(ds, (480,480,1,1))


savecube(ds, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/SIF_TROPOMI_gridded_10km_filtered_europe_2018-05-04_2020-12-26.zarr"; overwrite = true)




## Sentinel-3 LST

lst_files = glob("LST_*.nc", "/Net/Groups/BGI/scratch/jgens/Sentinel3/gridded/LST_1km/")

lst_files[1]

grid = Cube(open_dataset("/Net/Groups/BGI/scratch/jgens/Sentinel3/gridded/LST_1km/grid_h17-h20_v2-v5_1km.nc"))

lat = grid.data[:,:,1]'


lon = grid.data[:,:,2]'


# time dimension

time = lst_files

time = DateTime.([time[i][end-12:end-3] for i in eachindex(time)])

# variable dimension

first_lst = Cube(open_dataset(lst_files[1]))

Variable = collect(lookup(first_lst, :Variable))

push!(Variable, "lon", "lat")

# YAXArrays dimension

data_empty = Array{Float32, 4}(undef, 4800,4800,length(time),4)

axlist = (
    Dim{:X}(1:4800),
    Dim{:Y}(1:4800),
    Dim{:time}(time),
    Dim{:Variable}(Variable))




for i in eachindex(lst_files)

    tmp_cube = Cube(open_dataset(lst_files[i]))
    new_cube = collect(tmp_cube.data)

    for d in 1:2
        new_cube[:,:,d] = reverse(new_cube[:,:,d], dims=2)'
    end

    data_empty[:,:,i,1:2] = new_cube

    data_empty[:,:,i,3] = lon

    data_empty[:,:,i,4] = lat

end


props = Dict(
    "X" => "lon",
    "Y" => "lat",
    "time" => "days",
    "lon" => "longitude",
    "lat" => "latitude",
    "gridding" => "1km",
    "proj" => "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
    "Author" => "MPI-BGC, Daniel Pabon (dpabon@bgc-jena.mpg.de)"
)

ds = YAXArray(axlist, data_empty, props)

ds = setchunks(ds, (4800,4800,1,1))


savecube(ds, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/LST_Sentinel3_gridded_1km_europe_2018-01-04_2020-12-29.zarr"; overwrite = true)


## Sentinel-3 OLCI

OLCI_files = glob("OLCI_*.nc", "/Net/Groups/BGI/scratch/jgens/Sentinel3/gridded/OLCI_1km/")

OLCI_files[1]

grid = Cube(open_dataset("/Net/Groups/BGI/scratch/jgens/Sentinel3/gridded/OLCI_1km/grid_h17-h20_v2-v5_1km.nc"))

lat = grid.data[:,:,1]'


lon = grid.data[:,:,2]'


# time dimension

time = OLCI_files

time = DateTime.([time[i][end-12:end-3] for i in eachindex(time)])

# variable dimension

first_OLCI = Cube(open_dataset(OLCI_files[1]))

Variable = collect(lookup(first_OLCI, :Variable))

push!(Variable, "lon", "lat")

# YAXArrays dimension

data_empty = Array{Float32, 4}(undef, 4800,4800,length(time),10)

axlist = (
    Dim{:X}(1:4800),
    Dim{:Y}(1:4800),
    Dim{:time}(time),
    Dim{:Variable}(Variable))




for i in eachindex(OLCI_files)

    tmp_cube = Cube(open_dataset(OLCI_files[i]))
    new_cube = collect(tmp_cube.data)

    for d in 1:8
        new_cube[:,:,d] = reverse(new_cube[:,:,d], dims=2)'
    end

    data_empty[:,:,i,1:8] = new_cube

    data_empty[:,:,i,9] = lon

    data_empty[:,:,i,10] = lat

end


props = Dict(
    "X" => "lon",
    "Y" => "lat",
    "time" => "days",
    "lon" => "longitude",
    "lat" => "latitude",
    "gridding" => "1km",
    "proj" => "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
    "Author" => "MPI-BGC, Daniel Pabon (dpabon@bgc-jena.mpg.de)"
)

ds = YAXArray(axlist, data_empty, props)

ds = setchunks(ds, (4800,4800,1,1))


savecube(ds, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/OLCI_Sentinel3_gridded_1km_europe_2018-01-04_2020-12-29.zarr"; overwrite = true)



### MODIS MC43A4-1km ###


MC43A4_files = glob("MCD43A4_*.nc", "/Net/Groups/BGI/scratch/jgens/Modis/gridded/MCD43A4_1km/")

open_mfdataset(MC43A4_files)

YAXArrays.Datasets.open_mfdataset(MC43A4_files)

YAXArrays.Datasets.open_mfdataset("/Net/Groups/BGI/scratch/jgens/Modis/gridded/MCD43A4_1km/MCD*.nc")

MC43A4_files[1]

# Extracting latitude and longitude vectors

grid = Cube(open_dataset("/Net/Groups/BGI/scratch/jgens/Modis/gridded/MCD43A4_1km/grid_h17-h20_v2-v5_1km.nc"))

grid.data[:,:,1]


lat = grid.data[:,:,1]'

lon = grid.data[:,:,2]'

# time dimension

time = MC43A4_files

time = DateTime.([time[i][end-12:end-3] for i in eachindex(time)])

# variable dimension

first_MC43A4 = Cube(open_dataset(MC43A4_files[1]))

Variable = collect(lookup(first_MC43A4, :bands))

Variable = push!(string.(Variable), "lon", "lat")

# YAXArrays dimension


axlist = (
Dim{:X}(1:4800),
Dim{:Y}(1:4800),
Dim{:time}(time),
Dim{:Variable}(Variable))

z1 = zcreate(Float32, 4800,4800,976,9,path = "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/MC43A4_gridded_1km_europe_2018-05-01_2020-12-31.zarr",chunks=(4800, 4800,1,1))

z1


@showprogress for i in eachindex(MC43A4_files)
    tmp_cube = Cube(open_dataset(MC43A4_files[i]))
    new_cube = collect(tmp_cube.data)

    for d in 1:7
        new_cube[:,:,d] = reverse(new_cube[:,:,d], dims=2)'
    end

    z1[:,:,i,1:7] = new_cube

    z1[:,:,i,8] = lon

    z1[:,:,i,9] = lat

end


props = Dict(
    "X" => "lon",
    "Y" => "lat",
    "time" => "days",
    "lon" => "longitude",
    "lat" => "latitude",
    "gridding" => "10km",
    "proj" => "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
    "Author" => "MPI-BGC, Daniel Pabon (dpabon@bgc-jena.mpg.de)"
)

ds = YAXArray(axlist, z1, props)

ds = setchunks(ds, (4800,4800,1,1))


savecube(ds, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/MC43A4_gridded_1km_europe_2018-05-01_2020-12-31_YAXA.zarr")

#############################################


# MOD11A1

MOD11A1_files = glob("MOD11A1_*.nc", "/Net/Groups/BGI/scratch/jgens/Modis/gridded/MOD11A1_1km/")


MOD11A1_files[1]

# Extracting latitude and longitude vectors

grid = Cube(open_dataset("/Net/Groups/BGI/scratch/jgens/Modis/gridded/MOD11A1_1km/grid_h17-h20_v2-v5_1km.nc"))

grid.data[:,:,1]


lat = grid.data[:,:,1]'

lon = grid.data[:,:,2]'

# time dimension

time = MOD11A1_files

time = DateTime.([time[i][end-21:end-12] for i in eachindex(time)])

# variable dimension

first_MOD11A1 = Cube(open_dataset(MOD11A1_files[1]))

Variable = collect(lookup(first_MOD11A1, :Variable))

Variable = push!(string.(Variable), "lon", "lat")

# YAXArrays dimension


axlist = (
Dim{:X}(1:4800),
Dim{:Y}(1:4800),
Dim{:time}(time),
Dim{:Variable}(Variable))

z1 = zcreate(Float32, 4800,4800,length(time),4,path = "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/MOD11A1_gridded_1km_europe_2018-05-01_2020-12-31.zarr",chunks=(4800, 4800,1,1))

z1


#data_empty = Array{Float32, 4}(undef, 4800,4800,length(time),4)


@showprogress for i in eachindex(MOD11A1_files)
    tmp_cube = Cube(open_dataset(MOD11A1_files[i]))
    new_cube = collect(tmp_cube.data)

    for d in 1:2
        new_cube[:,:,d] = reverse(new_cube[:,:,d], dims=2)'
    end

    z1[:,:,i,1:2] = new_cube

    z1[:,:,i,3] = lon

    z1[:,:,i,4] = lat

end


props = Dict(
    "X" => "lon",
    "Y" => "lat",
    "time" => "days",
    "lon" => "longitude",
    "lat" => "latitude",
    "gridding" => "10km",
    "proj" => "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
    "Author" => "MPI-BGC, Daniel Pabon (dpabon@bgc-jena.mpg.de)"
)

ds = YAXArray(axlist, z1, props)

ds = setchunks(ds, (4800,4800,1,1))


savecube(ds, "/Net/Groups/BGI/work_3/OEMC/oemc_sif/data/MOD11A1_gridded_1km_europe_2018-05-01_2020-12-31_YAXA.zarr")

