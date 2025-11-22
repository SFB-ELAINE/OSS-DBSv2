import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("Results/stimulation_in_time.csv")
plt.title("Should be zero")
plt.plot(data["time"] / 1e-3, data["E1C1"], label="Applied")
plt.plot(data["time"] / 1e-3, data["E1C1_free"], label="Free")
plt.legend()
plt.show()
plt.title("Should be the applied signal")
plt.plot(data["time"] / 1e-3, data["E1C2"], label="Applied")
plt.plot(data["time"] / 1e-3, data["E1C2_free"], label="Free")
plt.legend()
plt.show()
