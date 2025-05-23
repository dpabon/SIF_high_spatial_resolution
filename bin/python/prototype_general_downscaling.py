import xarray as xr
import numpy as np
import xesmf as xe
from scipy.optimize import minimize
import inspect
import os # For path operations

def downscaling(dataset_low_res: xr.Dataset,
                dataset_high_res: xr.Dataset,
                model,  # Model function, e.g., model(params_array, pred1_array, pred2_array, ...)
                cost_f, # Cost function, e.g., cost_f(params, pred1_data, ..., predN_data, obs_data)
                param_ini: np.ndarray, # Initial parameters for optimization, 1D array
                param_bounds: np.ndarray, # Bounds for parameters, shape (2, num_params)
                min_obs: int = 25,
                winsize_opt: int = 5,
                save_inter_folder: str = None,
                method_regriding: str = 'bilinear',
                optimizer_method: str = 'L-BFGS-B', 
                optimizer_options: dict = None):    
    """
    Performs downscaling by optimizing parameters of a given model using a moving window approach.

    Args:
        dataset_low_res: xarray.Dataset containing the low-resolution observed variable (target).
                         Must contain exactly one data variable.
        dataset_high_res: xarray.Dataset containing high-resolution predictor variables.
        model: Callable model function. Expected signature: model(params, pred_var1, ..., pred_varN),
               where params is an array of optimizable parameters, and pred_varX are predictor data arrays.
        cost_f: Callable cost function to be minimized. Expected signature:
                cost_f(optim_params, pred1_data_flat, ..., predN_data_flat, obs_data_flat).
                It should return a scalar value.
        param_ini: Initial guess for the parameters to be optimized (1D NumPy array).
        param_bounds: Bounds for the parameters (NumPy array of shape (2, num_parameters)),
                      where param_bounds[0,:] are lower bounds and param_bounds[1,:] are upper bounds.
        min_obs: Minimum number of valid observations in a window to perform optimization.
        winsize_opt: Size of the moving window for optimization (applied to lat and lon).
        save_inter_folder: Optional. Path to a folder to save intermediate regridded dataset.
        method_regriding: Method for regridding (e.g., 'bilinear', 'conservative').
        optimizer_method: Method to be used by scipy.optimize.minimize (e.g., 'L-BFGS-B', 'SLSQP').
        optimizer_options: Dictionary of options for the optimizer (e.g., {'maxiter': 1000, 'ftol': 1e-7}).
                           If None, defaults to {'maxiter': 1000, 'ftol': 1e-7, 'gtol': 1e-5} for L-BFGS-B.

    Returns:
        xarray.DataArray: Cube of optimized parameters.
    """

    # --- Input Validations ---
    # Check 1: dataset_low_res has only one variable
    low_res_vars = list(dataset_low_res.data_vars.keys())
    if len(low_res_vars) != 1:
        raise ValueError('dataset_low_res must contain exactly one data variable.')
    dataset_low_res_var_name = low_res_vars[0]

    # Check 2: Number of predictor variables in model signature matches dataset_high_res
    try:
        model_signature = inspect.signature(model)
        num_model_data_args = len(model_signature.parameters) - 1 
    except TypeError:
        raise TypeError("The 'model' argument must be a callable function or method.")
        
    high_res_vars = list(dataset_high_res.data_vars.keys())
    if num_model_data_args < 0 : 
        raise ValueError(
            f"The model function '{model.__name__}' signature is not as expected. "
            "It should take optimizable parameters as its first argument."
        )
    if num_model_data_args != len(high_res_vars):
        raise ValueError(
            f"The model function '{model.__name__}' expects {num_model_data_args} data (predictor) variables "
            f"(after the initial parameter argument), but dataset_high_res provides {len(high_res_vars)} variables. "
            "These numbers must match. Adjust the model definition or dataset_high_res."
        )

    # Check 3: param_ini and param_bounds dimensions
    if not isinstance(param_ini, np.ndarray) or param_ini.ndim != 1:
        raise ValueError("param_ini must be a 1D NumPy array.")
    num_optim_params = len(param_ini)
    
    if not isinstance(param_bounds, np.ndarray) or param_bounds.shape != (2, num_optim_params):
        raise ValueError(
            f"param_bounds must be a NumPy array of shape (2, {num_optim_params}), "
            f"but got shape {param_bounds.shape}."
        )

    # Set default optimizer options if None
    if optimizer_options is None:
        optimizer_options = {'maxiter': 1000, 'ftol': 1e-7, 'gtol': 1e-5}


    # --- Regridding high-res data to low-res grid ---
    regridder = xe.Regridder(dataset_high_res, dataset_low_res, method=method_regriding)
    dataset_low_pred = regridder(dataset_high_res) 

    if save_inter_folder is not None:
        if not isinstance(save_inter_folder, str):
            raise TypeError("save_inter_folder must be a string path if provided.")
        os.makedirs(save_inter_folder, exist_ok=True) 
        
        try:
            dataset_low_pred_chunked = dataset_low_pred.chunk(dataset_low_res.chunks)
        except (ValueError, TypeError) as e:
            print(f"Warning: Could not chunk dataset_low_pred using dataset_low_res.chunks directly: {e}. "
                  "Using auto chunking for lat/lon/time dimensions if present.")
            auto_chunks = {dim: 'auto' for dim in dataset_low_pred.dims if dim in ['lat', 'lon', 'time']}
            if not auto_chunks: 
                 auto_chunks = {dim: 'auto' for dim in dataset_low_pred.dims}
            if not auto_chunks: 
                dataset_low_pred_chunked = dataset_low_pred 
            else:
                dataset_low_pred_chunked = dataset_low_pred.chunk(auto_chunks)

        zarr_filename = 'intermediate_results_low_resolution_regridded_dataset.zarr'
        zarr_path = os.path.join(save_inter_folder, zarr_filename)
        dataset_low_pred_chunked.to_zarr(zarr_path, mode='w') 
        dataset_low_pred = xr.open_zarr(zarr_path, chunks={}) 

    # --- Optimizing parameters using a moving window ---
    window_size_lat = winsize_opt
    window_size_lon = winsize_opt

    # Generalized optimization function for a single window
    def optimize_params_window_generalized(obs_cube_window, *pred_cubes_window_tuple,
                                           param_ini_optim, param_bounds_optim, 
                                           min_obs_val, cost_func_optim,
                                           current_optimizer_method, current_optimizer_options): # Added optimizer params
        """
        Wrapper function to optimize parameters for a single window.
        Designed for use with xr.apply_ufunc.
        """
        all_cubes_in_window = [obs_cube_window] + list(pred_cubes_window_tuple)
        
        mask = np.ones_like(all_cubes_in_window[0], dtype=bool) 
        for cube_data in all_cubes_in_window:
            mask &= ~np.isnan(cube_data)
        
        n_valid = np.sum(mask)

        if n_valid < min_obs_val:
            return np.full(len(param_ini_optim), np.nan, dtype=np.float32)

        obs_flat = obs_cube_window[mask]
        predictors_flat = [pred_cube[mask] for pred_cube in pred_cubes_window_tuple]
        
        bounds_scipy = list(zip(param_bounds_optim[0, :], param_bounds_optim[1, :]))
        cost_args = tuple(predictors_flat) + (obs_flat,)
        
        try:
            result = minimize(cost_func_optim,
                              x0=param_ini_optim,
                              args=cost_args,
                              method=current_optimizer_method, # Use passed method
                              bounds=bounds_scipy,
                              options=current_optimizer_options) # Use passed options

            if result.success:
                optimized_params = np.clip(result.x, param_bounds_optim[0, :], param_bounds_optim[1, :])
            else:
                optimized_params = np.full(len(param_ini_optim), np.nan)
        except Exception as e:
            optimized_params = np.full(len(param_ini_optim), np.nan)
            
        return optimized_params.astype(np.float32) 

    # --- Prepare inputs for xr.apply_ufunc ---
    rolling_obs_var = dataset_low_res[dataset_low_res_var_name].rolling(
        lat=window_size_lat, lon=window_size_lon, center=True, min_periods=1 
    ).construct(lat='lat_roll', lon='lon_roll')

    input_arrays_for_ufunc = [rolling_obs_var]
    
    dataset_low_pred_var_names = list(dataset_low_pred.data_vars.keys()) 
    for var_name in dataset_low_pred_var_names:
        rolling_pred_var = dataset_low_pred[var_name].rolling(
            lat=window_size_lat, lon=window_size_lon, center=True, min_periods=1
        ).construct(lat='lat_roll', lon='lon_roll')
        input_arrays_for_ufunc.append(rolling_pred_var)

    input_core_dims_for_ufunc = [['lat_roll', 'lon_roll']] * len(input_arrays_for_ufunc)
    
    # --- Perform the optimization using apply_ufunc ---
    parameters_cube = xr.apply_ufunc(
        optimize_params_window_generalized, 
        *input_arrays_for_ufunc,          
        kwargs={                          
            'param_ini_optim': param_ini,
            'param_bounds_optim': param_bounds,
            'min_obs_val': min_obs,
            'cost_func_optim': cost_f,
            'current_optimizer_method': optimizer_method, # Pass new optimizer method
            'current_optimizer_options': optimizer_options # Pass new optimizer options
        },
        input_core_dims=input_core_dims_for_ufunc, 
        output_core_dims=[['parameters']],       
        dask_gufunc_kwargs={'output_sizes': {'parameters': num_optim_params}}, 
        dask='parallelized',                
        output_dtypes=[np.float32],         
        vectorize=True,                     
        exclude_dims=set(('lat_roll', 'lon_roll')), 
    ).rename('optimized_parameters')
    
    return parameters_cube

