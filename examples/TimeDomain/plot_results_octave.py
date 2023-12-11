import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("ResultsOctave/stimulation_in_time.csv")


fig, ax1 = plt.subplots()
plt.title("Should be zero")
color = "tab:red"

ax1.set_xlabel("Time / ms")
ax1.plot(data["time"] / 1e-3, data["E1C1"], label="Applied", color=color)
ax1.tick_params(axis="y", labelcolor=color)

ax2 = ax1.twinx()
color = "tab:blue"
ax2.plot(data["time"] / 1e-3, data["E1C1_free"], label="Free", color=color)
ax2.tick_params(axis="y", labelcolor=color)

fig.tight_layout()
plt.legend()
plt.show()

fig, ax1 = plt.subplots()
plt.title("Should be the applied signal")
color = "tab:red"

ax1.set_xlabel("Time / ms")
ax1.plot(data["time"] / 1e-3, data["E1C2"], label="Applied", color=color)
plt.legend()
ax1.tick_params(axis="y", labelcolor=color)

ax2 = ax1.twinx()
color = "tab:blue"
ax2.plot(data["time"] / 1e-3, data["E1C2_free"], label="Free", color=color)
ax2.tick_params(axis="y", labelcolor=color)

plt.legend()
fig.tight_layout()
plt.show()
