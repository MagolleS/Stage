from equation1_euler import *
from equation2_euler import *
from equation1_odeint import *
from equation2_odeint import *
from equation2_euler_2d import *
from fonctions import *
import random
import pickle

################Barres 3D ################
# nb_grad_x=100
# S_0 = np.array([random.random()*10 for x in range(nb_grad_x)])
# N1_0 = np.array([random.random()*0 for x in range(nb_grad_x)])
# N2_0 = np.array([random.random()*10 for x in range(nb_grad_x)])
# R_S = 3
# D_N = 2.5/10000
# S, N1, N2, temps = resolve2_euler((S_0, N1_0, N2_0), 200, 50, 100, D_N=D_N, R_S=R_S, proba_ressource=0.005, proba_N1=0.001, R_N1=100, arrondir_cellules=True)
# N1_i, N2_i = integrate(N1, 1/100), integrate(N2, 1/100)
# N1_m, N2_m = moyenne(N1_i[150*1000:]), moyenne(N2_i[150*1000:])
# print("N1_m = {}, N2_m = {}".format(N1_m, N2_m))



# nb_grad_x = 100
# list_R_S = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
# liste_D_N = [25/10000, 10/10000,  5/10000, 2.5/10000]


# l_N1_top = []
# l_N2_top = []

# for D_N in liste_D_N:
# 	liste_int_1 = []
# 	liste_int_2 = []
# 	for R_S in list_R_S:
# 		print("D_N = {}, R_S = {}".format(D_N, R_S))
# 		S_0 = np.array([random.random()*10 for x in range(nb_grad_x)])
# 		N1_0 = np.array([random.random()*0 for x in range(nb_grad_x)])
# 		N2_0 = np.array([random.random()*10 for x in range(nb_grad_x)])

# 		S, N1, N2, temps = resolve2_euler((S_0, N1_0, N2_0), 200, 50, 100, D_N=D_N, R_S=R_S, proba_ressource=0.005, proba_N1=0.001, R_N1=100, arrondir_cellules=True)
# 		N1_i, N2_i = integrate(N1, 1/100), integrate(N2, 1/100)
# 		N1_m, N2_m = moyenne(N1_i[150*1000:]), moyenne(N2_i[150*1000:])


# 		liste_int_1.append(N1_m)
# 		liste_int_2.append(N2_m)

# 	l_N1_top = liste_int_1 + l_N1_top
# 	l_N2_top = liste_int_2 + l_N2_top

# with open("fig_3d.data", "rb") as f:
# 	l_N1_top, l_N2_top = pickle.load(f)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")
# # ax.set_ylim([2.5/10000, 5/10000])

# ax.set_yticks(range(len(liste_D_N)))
# ax.set_yticklabels(liste_D_N[::-1])
# ax.invert_yaxis()

# ax.set_xlabel("Ressource influx")
# ax.set_ylabel("Cell diffusion")
# ax.set_zlabel("Population size")
# x_, y_  = np.meshgrid(list_R_S, range(len(liste_D_N)))
# x, y = x_.ravel(), y_.ravel()

# line_red, line_blue = Patch(facecolor='red', edgecolor="black"), Patch(facecolor='blue', edgecolor="black")
# plt.legend([line_blue, line_red], ["N1", "N2"])
# witdh, depth = 0.15, 0.25
# ax.bar3d(x, y, l_N1_top, witdh, depth, l_N2_top, color="red", shade=False, edgecolor="black")
# ax.bar3d(x, y, [0]*len(l_N2_top), witdh, depth, l_N1_top, color="blue", shade=False, edgecolor="black")
# ax.set_zlim([0,200])

# plt.show()
################2D ################
n = 100*100
S_0  = np.array([random.random()*10 for i in range(n)])
N1_0  = np.array([0*random.random() for i in range(n)])
N2_0  = np.array([random.random()*10 for i in range(n)])

S, N1, N2, temps = resolve2_euler_2d((S_0, N1_0, N2_0), 10, 3, 6, R_S = 15, proba_ressource = 0.005, R_N1 = 100, proba_N1=0.01, arrondir_cellules = True)
afficher2D(S, N1, N2, temps, 100,100,show_integrate=True)


################Comparaison 3 ################
# nb_grad_x = 100
# dx = 1/100

# S_0 = np.array([random.random()*10 for x in range(nb_grad_x)])
# N1_0 = np.array([random.random()*0 for x in range(nb_grad_x)])
# N2_0 = np.array([random.random()*10 for x in range(nb_grad_x)])

# S_sa, N1_sa, N2_sa, temps1 = resolve2_euler((S_0, N1_0, N2_0), 150, 50, 100, R_S=1.5, proba_ressource= 0.005, R_N1 = 15, proba_N1=0.001, arrondir_cellules=False )
# S_a, N1_a, N2_a, temps2 = resolve2_euler((S_0, N1_0, N2_0), 150, 50, 100, R_S=1.5, proba_ressource= 0.005, R_N1 = 15, proba_N1=0.001, arrondir_cellules=True)

# S_sa, N1_sa, N2_sa = integrate(S_sa, dx), integrate(N1_sa, dx), integrate(N2_sa, dx)
# S_a, N1_a, N2_a = integrate(S_a, dx), integrate(N1_a, dx), integrate(N2_a, dx)

# fig, ax = plt.subplots()
# ax.set_ylabel("Population size")
# ax.set_xlabel("t")
# ax.tick_params(axis="y")

# ax_S = ax.twinx()
# ax_S.set_ylabel("Ressource size")
# ax_S.tick_params(axis="y", labelcolor="green")


# ax_S.plot(temps1, S_sa, label='S', color="green")
# ax.plot(temps1, N1_sa, label='N1', color = "blue")
# ax.plot(temps1, N2_sa, label='N2', color = "red")


# ax.legend(loc = "upper right")
# ax_S.legend(loc = "center right")
# plt.title("Sans arrondir")



# fig, ax = plt.subplots()
# ax.set_ylabel("Population size")
# ax.set_xlabel("t")
# ax.tick_params(axis="y")

# ax_S = ax.twinx()
# ax_S.set_ylabel("Ressource size")
# ax_S.tick_params(axis="y", labelcolor="green")



# ax_S.plot(temps2, S_a, label='S', color="green")
# ax.plot(temps2, N1_a, label='N1', color = "blue")
# ax.plot(temps2, N2_a, label='N2', color = "red")



# ax.legend(loc = "upper right")
# ax_S.legend(loc = "center right")
# plt.title("En arrondissant")

# plt.show()




################Comparaison 2 ################
# nb_grad_x = 100
# dx = 1/100

# S_0 = np.array([random.random()*1 for x in range(nb_grad_x)])
# N1_0 = np.array([random.random()*0 for x in range(nb_grad_x)])
# N2_0 = np.array([random.random()*1 for x in range(nb_grad_x)])

# S_sa, N1_sa, N2_sa, temps1 = resolve2_euler((S_0, N1_0, N2_0), 150, 50, 100, R_S=1.5, proba_ressource= 0.005, R_N1 = 30, proba_N1=0.001, arrondir_cellules=False )
# S_a, N1_a, N2_a, temps2 = resolve2_euler((S_0, N1_0, N2_0), 150, 50, 100, R_S=1.5, proba_ressource= 0.005, R_N1 = 30, proba_N1=0.001, arrondir_cellules=True)

# S_sa, N1_sa, N2_sa = integrate(S_sa, dx), integrate(N1_sa, dx), integrate(N2_sa, dx)
# S_a, N1_a, N2_a = integrate(S_a, dx), integrate(N1_a, dx), integrate(N2_a, dx)

# fig, ax = plt.subplots()
# ax.set_ylabel("Population size")
# ax.set_xlabel("t")
# ax.tick_params(axis="y")

# ax_S = ax.twinx()
# ax_S.set_ylabel("Ressource size")
# ax_S.tick_params(axis="y", labelcolor="green")
# ax_S.set_ylim([0,1.6])


# ax_S.plot(temps1, S_sa, label='S', color="green")
# ax.plot(temps1, N1_sa, label='N1', color = "blue")
# ax.plot(temps1, N2_sa, label='N2', color = "red")


# ax.legend()
# ax_S.legend(loc = "center right")
# plt.title("Sans arrondir")



# fig, ax = plt.subplots()
# ax.set_ylabel("Population size")
# ax.set_xlabel("t")
# ax.tick_params(axis="y")

# ax_S = ax.twinx()
# ax_S.set_ylabel("Ressource size")
# ax_S.tick_params(axis="y", labelcolor="green")
# ax_S.set_ylim([0,1.6])



# ax_S.plot(temps2, S_a, label='S', color="green")
# ax.plot(temps2, N1_a, label='N1', color = "blue")
# ax.plot(temps2, N2_a, label='N2', color = "red")



# ax.legend(loc = "upper right")
# ax_S.legend(loc = "center right")
# plt.title("En arrondissant")

# plt.show()



################Comparaison 1 ################
# S_0, N1_0, N2_0  = 5,10,5
# nb_grad_x = 100

# S_0_ = np.array([random.random()*(S_0*2) for x in range(nb_grad_x)])
# N1_0_ = np.array([random.random()*(N1_0*2) for x in range(nb_grad_x)])
# N2_0_ = np.array([random.random()*(N2_0*2) for x in range(nb_grad_x)])

# print("Solving 1.5s sol1...")
# S_sd_1, N1_sd_1, N2_sd_1, temps1 = resolve1_odeint((S_0, N1_0, N2_0), 1.5)

# print("Solving 50s sol1...")
# S_sd_50, N1_sd_50, N2_sd_50, temps2 = resolve1_odeint((S_0, N1_0, N2_0), 50)


# print("Solving 1.5s sol2...")
# S_d_1, N1_d_1, N2_d_1, temps3 = resolve2_odeint((S_0_, N1_0_, N2_0_), 1.5)
# S_d_1, N1_d_1, N2_d_1 = integrate(S_d_1, 1/100), integrate(N1_d_1, 1/100), integrate(N2_d_1, 1/100)


# print("Solving 50s sol2...")
# S_d_50, N1_d_50, N2_d_50, temps4 = resolve2_odeint((S_0_, N1_0_, N2_0_), 50)
# S_d_50, N1_d_50, N2_d_50 = integrate(S_d_50, 1/100), integrate(N1_d_50, 1/100), integrate(N2_d_50, 1/100)


# print(len(temps1))
# print(len(temps3))
# print()

# print(len(temps2))
# print(len(temps4))


# plt.figure()
# fig, ax = plt.subplots()
# ax.set_ylabel("Population size")
# ax.set_xlabel("t")
# ax.tick_params(axis="y")

# ax_S = ax.twinx()
# ax_S.set_ylabel("Ressource size")
# ax_S.tick_params(axis="y", labelcolor="green")
# ax_S.set_ylim([0,1])


# ax_S.plot(temps1, S_d_1, label='S', color="green")
# ax_S.plot(temps1, S_sd_1, '--', label='S uniforme', color="green")

# ax.plot(temps1, N1_d_1, label='N1', color = "blue")
# ax.plot(temps1, N1_sd_1, '--', label='N1 uniforme', color = "blue")

# ax.plot(temps1, N2_d_1, label='N2', color = "red")
# ax.plot(temps1, N2_sd_1, '--', label='N2 uniforme', color = "red")



# ax.legend()
# ax_S.legend(loc = "center right")



# plt.figure()
# fig, ax = plt.subplots()
# ax.set_ylabel("Population size")
# ax.set_xlabel("t")
# ax.tick_params(axis="y")

# ax_S = ax.twinx()
# ax_S.set_ylabel("Ressource size")
# ax_S.tick_params(axis="y", labelcolor="green")
# ax_S.set_ylim([0,1])


# ax_S.plot(temps2, S_d_50, label='S', color="green")
# ax_S.plot(temps2, S_sd_50, '--', label='S uniforme', color="green")

# ax.plot(temps2, N1_d_50, label='N1', color = "blue")
# ax.plot(temps2, N1_sd_50, '--', label='N1 uniforme', color = "blue")

# ax.plot(temps2, N2_d_50, label='N2', color = "red")
# ax.plot(temps2, N2_sd_50, '--', label='N2 uniforme', color = "red")



# ax.legend(loc = "upper right")
# ax_S.legend(loc = "center right")

# plt.show()