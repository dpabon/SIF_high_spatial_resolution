using Pkg
Pkg.activate("/Net/Groups/BGI/work_3/OEMC/oemc_sif/bin/")
using YAXArrays, NetCDF, Zarr



PathLoad = "/Net/Groups/BGI/scratch/jgens/Sen4GPP/Fluxcom_SIF/"

cube_loc = "/Net/Groups/BGI/scratch/jgens/FLUXCOM_TROPOSIF"

START_YEAR = 2018
END_YEAR   = 2021

PRODUCT = "TROPOSIF"

sites = ["BE-Bra", "BE-Dor", "BE-Lcr", "BE-Lon", "BE-Maa", "BE-Vie",
         "CH-Aws", "CH-Cha", "CH-Dav", "CH-Fru", "CH-Lae", "CH-Oe2",
         "CZ-BK1", "CZ-KrP", "CZ-Lnz", "CZ-RAJ", "CZ-Stn", "CZ-wet",
         "DE-Akm", "DE-Geb", "DE-Gri", "DE-Hai", "DE-HoH", "DE-Hzd",
         "DE-Kli", "DE-Obe", "DE-RuR", "DE-RuS", "DE-RuW", "DE-Tha",
         "DK-Gds", "DK-Sor", "ES-Abr", "ES-Agu", "ES-Cnd", "ES-LJu",
         "ES-LM1", "ES-LM2", "FI-Hyy", "FI-Ken", "FI-Let", "FI-Qvd",
         "FI-Sii", "FI-Var", "FR-Aur", "FR-Bil", "FR-FBn", "FR-Fon",
         "FR-Gri", "FR-Hes", "FR-LGt", "FR-Lam", "FR-Tou", "GF-Guy",
         "GL-Dsk", "IE-Cra", "IL-Yat", "IT-BCi", "IT-BFt", "IT-Cp2",
         "IT-Lav", "IT-Lsn", "IT-MBo", "IT-Ren", "IT-SR2", "IT-Tor",
         "RU-Fy2", "RU-Fyo", "SE-Deg", "SE-Htm", "SE-Nor", "SE-Ros", "SE-Svb"]

sites = ["BE-Bra"]

test_cube = open_dataset("/Net/Groups/BGI/scratch/jgens/Sen4GPP/Fluxcom_SIF/BE-Bra/BE-Bra_TROPOSIF_2018.nc")
test_cube.time

length(test_cube.time) 

open_dataset("/Net/Groups/BGI/scratch/jgens/TROPOSIF/gridded/global/TROPOSIF_005_2018.nc")

