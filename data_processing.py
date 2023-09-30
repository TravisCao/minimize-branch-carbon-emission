from scipy.io import loadmat, savemat
import numpy as np

# solar
solar_data = loadmat("data/texas_2020_solar.mat")["solar_MW"]
solar_pattern = np.mean(solar_data.reshape(366, 24, -1), axis=(0, 2))
# save solar as mat file
solar_dict = {"solar": solar_pattern}
savemat("./solar.mat", solar_dict)

wind_data = loadmat("data/texas_2020_wind.mat")["wind_MW"]
wind_pattern = np.mean(wind_data.reshape(366, 24, -1), axis=(0, 2))
# save wind as mat file
wind_dict = {"wind": wind_pattern}
savemat("./wind.mat", wind_dict)

hydro_data = loadmat("data/texas_2020_hydro.mat")["hydro_MW"]
hydro_pattern = np.mean(hydro_data.reshape(366, 24, -1), axis=(0, 2))
# save hydro as mat file
hydro_dict = {"hydro": hydro_pattern}
savemat("./hydro.mat", hydro_dict)


# plot wind_pattern
import matplotlib.pyplot as plt

plt.plot(wind_pattern)
plt.show()

plt.plot(solar_pattern)
plt.show()

plt.plot(hydro_pattern)
plt.show()
