import pandas as pd
from matplotlib import pyplot as plt 

output_P1 = pd.read_csv("output_P1.csv")

fig, axs = plt.subplots(2, 2)

axs[0, 0].plot(output_P1.time,output_P1.w_expl)
axs[0, 0].plot(output_P1.time,output_P1.exact_solution)
axs[0, 0].set_title("Euler forward method")
axs[0, 0].set_xlabel("t")
axs[0, 0].set_ylabel("y")
axs[0, 0].legend(["Euler forward","Exact solution"])

axs[1, 0].plot(output_P1.time,output_P1.w_impl)
axs[1, 0].plot(output_P1.time,output_P1.exact_solution)
axs[1, 0].set_title("Euler backward method")
axs[1, 0].set_xlabel("t")
axs[1, 0].set_ylabel("y")
axs[1, 0].legend(["Euler backward","Exact solution"])
axs[1, 0].sharex(axs[0, 0])
axs[1, 0].sharey(axs[0, 0])

axs[0, 1].plot(output_P1.time,output_P1.w_RK2)
axs[0, 1].plot(output_P1.time,output_P1.exact_solution)
axs[0, 1].set_title("Runge-Kutta with two stages (RK2) method")
axs[0, 1].set_xlabel("t")
axs[0, 1].set_ylabel("y")
axs[0, 1].legend(["Runge-Kutta","Exact solution"])
axs[0, 1].sharex(axs[0, 0])
axs[0, 1].sharey(axs[0, 0])

axs[1, 1].plot(output_P1.time,output_P1.w_RK2)
axs[1, 1].plot(output_P1.time,output_P1.exact_solution)
axs[1, 1].set_title("Crank-Nicolson (CN) method")
axs[1, 1].set_xlabel("t")
axs[1, 1].set_ylabel("y")
axs[1, 1].legend(["Crank-Nicolson","Exact solution"])
axs[1, 1].sharex(axs[0, 0])
axs[1, 1].sharey(axs[0, 0])

fig.tight_layout()
plt.show()



