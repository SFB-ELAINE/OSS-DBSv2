import matplotlib.pyplot as plt
import pandas as pd
data = pd.read_csv("Results/stimulation_in_time.csv")
plt.title("Should be zero")
plt.plot(data["time"], data["E1C1"])
plt.show()
plt.title("Should be the applied signal")
plt.plot(data["time"], data["E1C2"])
plt.show()
