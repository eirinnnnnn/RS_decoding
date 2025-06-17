import matplotlib.pyplot as plt
import pandas as pd
import io

# Load data from the user's CSV-like input
data = """
t0,t1,success_rate
32,0,1.000
30,1,1.000
28,2,1.000
26,3,1.000
24,4,1.000
22,5,1.000
20,6,0.985
18,7,0.989
16,8,0.996
14,9,0.999
12,10,1.000
"""

df = pd.read_csv(io.StringIO(data))

# Plotting
plt.figure(figsize=(8, 5))
plt.plot(df["t1"], df["success_rate"], marker='o')
plt.xlabel(r"$t_1$ (errors)")
plt.ylabel("Success Rate")
plt.title(r"False Acceptance Rate at $t_0 + 2t_1 = 32$")
plt.grid(True)
plt.ylim(0.98, 1.005)


plt.savefig("./rpt/fig/r32_rat.png")
