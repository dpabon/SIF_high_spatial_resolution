import xarray as xr
import numpy as np
import xesmf as xe
from scipy.optimize import minimize
import inspect
import os # For path operations


def downscaling(dataset_low_res: xr.Dataset, dataset_high_res: xr.Dataset, model, cost_f, param_ini = np.array([1.0, 2.0, -295.0, 10.0]), param_bounds = np.array([param_min, param_max]), min_obs = 25, winsize_opt = 5, save_inter_folder = None, method_regriding = 'bilinear'):


    # check that the number of variables in dataset_low_res is just one:

    if len(dataset_low_res.data_vars.keys()) != 1:
        raise IndexError('dataset_low_res contains more than one variable. Remove unnecessary variables and try again.')

    # check that the number of parameters in model is the same as the number of variables in the dataset

    if (len(inspect.signature(model).parameters.values()) - 1) != len(dataset_high_res.data_vars.keys()):
        raise IndexError('dataset_high_res contains different number of variables than necessary to run the model. Remove/add necessary variables, or modify the model function and try again.')
    
    # regridding to create low resolution dataset

    regridder = xe.Regridder(dataset_high_res, dataset_low_res, method=method_regriding)

    dataset_low_pred = regridder(dataset_high_res)

    if save_inter_folder != None:
        dataset_low_pred = dataset_low_pred.chunk(dataset_low_res.chunks)
        dataset_low_pred.to_zarr(save_inter_folder+'intermediate_results_low_resolution_regridded_dataset.zarr')

        dataset_low_pred = xr.open_dataset(save_inter_folder+'intermediate_results_low_resolution_regridded_dataset.zarr', chunks = {})

    # optimizing parameters using a moving window

    # Define window parameters
    window_size_lat = winsize_opt
    window_size_lon = winsize_opt

    if len(dataset_high_res.data_vars.keys()) == 2:


        def optimize_params_window(cube_low, cube_pred_low1, cube_pred_low2,
                            param_ini, param_bounds, min_obs):
            """
            Wrapper function to optimize SIF parameters for a single window.
            Designed to be used with apply_ufunc or iteration.
            Assumes input arrays (sif_w, etc.) are 2D numpy arrays (window data).
            """
            # Flatten the window arrays and filter out NaNs
            # Important: Filter consistently across all variables
            mask = ~np.isnan(cube_low) & ~np.isnan(cube_pred_low1) & ~np.isnan(cube_pred_low2)
            n_valid = np.sum(mask)

            if n_valid < min_obs:
                # Not enough valid observations in the window
                return np.full(len(param_ini), np.nan, dtype=np.float32)
            

            sif_obs_f = cube_low[mask]
            vi_f = cube_pred_low2[mask]
            #et_f = ndwi_w[mask] # Using NDWI as proxy for water stress/ET effect
            lst_f = cube_pred_low2[mask]

            # Define bounds for L-BFGS-B
            bounds_scipy = list(zip(*param_bounds)) # [(min1, max1), (min2, max2), ...]

        
            result = minimize(cost_function,
                                    x0=param_ini,
                                    args=(vi_f, lst_f, sif_obs_f),
                                    method='L-BFGS-B',
                                    bounds=bounds_scipy,
                                    options={'maxiter': 1000, 'ftol': 1e-7, 'gtol': 1e-5}) # Adjust options as needed

            optimized_params = np.array(np.clip(result.x, param_bounds[0], param_bounds[1]))
            return optimized_params
        
        parameters_cube = xr.apply_ufunc(
            optimize_params_window, # Function to apply
            # Input arrays:
            dataset_low_res[list(dataset_low_res.data_vars.keys())[0]].rolling(lat=window_size_lat, lon=window_size_lon, center=True).construct(lat = 'lat_roll', lon = 'lon_roll'),

            dataset_low_pred[list(dataset_low_pred.data_vars.keys())[0]].rolling(lat=window_size_lat, lon=window_size_lon, center=True).construct(lat = 'lat_roll', lon = 'lon_roll'),

            dataset_low_pred[list(dataset_low_pred.data_vars.keys())[1]].rolling(lat=window_size_lat, lon=window_size_lon, center=True).construct(lat = 'lat_roll', lon = 'lon_roll'),

            # Keyword arguments for the function:
            kwargs={'param_ini': param_ini, 'param_bounds': param_bounds, 'min_obs': min_obs_optim},
            # Define input core dimensions (the window dimensions):
            input_core_dims=[['lat_roll', 'lon_roll'], ['lat_roll', 'lon_roll'], ['lat_roll', 'lon_roll']],
            # Define output core dimensions (the parameters dimension):
            output_core_dims=[['parameters']],
            dask_gufunc_kwargs={'output_sizes': {'parameters': len(param_ini)}},
            # Specify the parameters dimension coordinates:
            dask='parallelized', # Enable Dask parallelization
            output_dtypes=[np.float64], # Specify output data type
            vectorize = True,
            # Add a new dimension 'parameters' to the output
            exclude_dims=set(('lat_roll', 'lon_roll')), # Exclude window dims from output
        ).rename('optimized_parameters')

    else:
        parameters_cube = np.nan

    return parameters_cube

