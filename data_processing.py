from scipy.io import loadmat, savemat
import numpy as np

# solar
solar_data = loadmat("data/texas_2020_solar.mat")["solar_MW"]
solar_pattern = np.mean(solar_data.reshape(366, 24, -1), axis=(0, 2))

wind_data = loadmat("data/texas_2020_wind.mat")["wind_MW"]
wind_pattern = np.mean(wind_data.reshape(366, 24, -1), axis=(0, 2))

hydro_data = loadmat("data/texas_2020_hydro.mat")["hydro_MW"]
hydro_pattern = np.mean(hydro_data.reshape(366, 24, -1), axis=(0, 2))

# scale data to around mean value
solar_pattern = solar_pattern / np.mean(solar_pattern)
wind_pattern = wind_pattern / np.mean(wind_pattern)
hydro_pattern = hydro_pattern / np.mean(hydro_pattern)

# 200 is set to be the PMAX of the solar and wind plant
solar_pattern = solar_pattern * 200
wind_pattern = wind_pattern * 200

# 1040 is the PMAX of the hydro plant
hydro_pattern = hydro_pattern * 600

# lambda function to convert numpy array and save mat file
savemat_helper = lambda arr, value: savemat(f"data/{value}.mat", {f"{value}": arr})

savemat_helper(solar_pattern, "solar")
savemat_helper(wind_pattern, "wind")
savemat_helper(hydro_pattern, "hydro")

# plot wind_pattern
import matplotlib.pyplot as plt

plt.plot(wind_pattern)
plt.show()

plt.plot(solar_pattern)
plt.show()

plt.plot(hydro_pattern)
plt.show()
