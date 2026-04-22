import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("Results_floating/stimulation_in_time.csv")
floatdata = pd.read_csv("Results_floating/floating_in_time.csv")
data_encap = pd.read_csv("Results_floating_encap/stimulation_in_time.csv")
floatdata_encap = pd.read_csv("Results_floating_encap/floating_in_time.csv")

plt.plot(1e6 * data["time"][0:400], data["E1C1"][0:400], label="Stimulation")
plt.plot(
    1e6 * floatdata["time"][0:400], floatdata["E1C2"][0:400], label="Passive contact"
)
plt.plot(
    1e6 * floatdata_encap["time"][0:400],
    floatdata_encap["E1C2"][0:400],
    label="Passive contact w/ encapsulation",
)
plt.legend()
plt.xlabel(r"Time / $\mu$s")
plt.ylabel("Voltage / V")
plt.savefig("Compare_contacts.pdf")
plt.show()
