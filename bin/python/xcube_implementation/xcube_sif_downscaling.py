#%%
from shapely.geometry import collection
import xcube
from xcube.core.store import new_data_store, get_data_store_params_schema

#%%
store_params = get_data_store_params_schema("stac")
store_params
# %%
store = new_data_store("stac", url="https://data-portal.s5p-pal.com/api/s5p-l3")

#%% 
store.get_data_ids()

# %%
store.describe_data('sif')
# %%
search_params = store.get_search_params_schema()
#%%
search_params

# %%
descriptors = list(
    store.search_data(
        time_range=["2024-06-18", "2024-06-20"]
    )
)
[d.to_dict() for d in descriptors]

# %%
descriptors
# %%
