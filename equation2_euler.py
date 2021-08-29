import numpy as np 
from fonctions import *
from math import floor, ceil
from matplotlib.patches import Patch
import pickle
import random

""" 
Résolution de l'équation 2  en 1D, avec la méthode d'euler pour la variable temporelle et odeint pour la variable spatiale.
1 : respiration
2 : fermentation
"""

x_0, x_f = 0,1

def resolve2_euler(initial_values, t_f, t_N1_0, t_N1_f, nu=15, d=10, c=1, D_N=5/10000, D_S=0.25/10000, R_S=0, R_N1=0, proba_ressource=0, proba_N1=0, arrondir_cellules=False, cste_fct=(32, 30, 2), dt = 1/1000, dx = 1/100):
	"""initial_values : valeurs des solutions au temps t=0 : (S_0, N1_0, N_20)
	t_f : temps d'arrêt de la simulation
	d,c, nu, D_N, D_S : paramètres du système d'équa diff
	proba_ressource : probabilité d'ajouter une quantité R_S de ressource à chaque pas de temps et chaque pas d'espace
	proba_N1 : probabilité d'ajouter une quantité R_N1 de respirateurs à chaque pas de temps et chaque pas d'espace (dans l'intervalle de temps [t_N1_0, t_N1_f])
	cste_fct : constantes qui définissent les fonctions (y_1, y_21, y_22)
	dt : pas de temps à utiliser
	"""
	S_0, N1_0, N2_0 = initial_values
	y_1, y_21, y_22 = cste_fct

	nb_grad_t = int(t_f/dt)
	nb_grad_x = int((x_f - x_0)/dx)
	dx2 = dx**2
	coef_prec = dt/dx2

	nu_l = np.array([nu]*nb_grad_x)
	d_l = np.array([d]*nb_grad_x)
	c_l = np.array([c]*nb_grad_x)

	S, N1, N2 = [S_0], [N1_0], [N2_0]
	S_t, N1_t, N2_t = S_0, N1_0, N2_0

	def maj(S):
		J1_S = S/(1+S)
		J2_S = 100*S/(100+S)
		J1_ATP = y_1*J1_S
		J2_ATP = y_21*J1_S + y_22*J2_S

		return J1_S, J2_S, J1_ATP, J2_ATP



	A = laplacienPeriodic(nb_grad_x)	
	M_inv_S = np.linalg.inv((np.identity(nb_grad_x)+D_S*coef_prec*A))
	M_inv_N = np.linalg.inv((np.identity(nb_grad_x)+D_N*coef_prec*A))

	temps = [0]
	t = 0
	while True:
		if abs(t-int(t))<10**(-3):
			print(t)

		J1_S, J2_S, J1_ATP, J2_ATP = maj(S_t)

		S_influx = np.array([0]*nb_grad_x)
		if R_S != 0:
			for i in range(nb_grad_x):
				if random.random() < proba_ressource:
					# print("ajout")
					S_influx[i] += R_S

		N1_influx = np.array([0]*nb_grad_x)
		if R_N1 != 0:
			if t_N1_0 <= t <= t_N1_f:
				for i in range(nb_grad_x):
					if random.random() < proba_N1:
						N1_influx[i] += R_N1

		N2_influx = np.array([0]*nb_grad_x)		
		
		S_t = M_inv_S.dot(S_t + dt*(nu_l - N1_t*J1_S - N2_t*J2_S)) + S_influx
		N1_t = M_inv_N.dot(N1_t + dt*N1_t*(c_l*J1_ATP - d_l)) + N1_influx
		N2_t = M_inv_N.dot(N2_t + dt*N2_t*(c_l*J2_ATP - d_l)) + N2_influx

		# if max(np.max(abs(S_t-S[-1])), np.max(abs(N1_t-N1[-1])), np.max(abs(N2_t-N2[-1]))) < seuil:
		# 	break

		if arrondir_cellules:
			for i in range(nb_grad_x):
				N1_i = N1_t[i]
				N2_i = N2_t[i]
				if random.random() < (N1_i-floor(N1_i)):
					N1_t[i] = ceil(N1_i)
				else:
					N1_t[i] = floor(N1_i)

				if random.random() < (N2_i-floor(N2_i)):
					N2_t[i] = ceil(N2_i)
				else:
					N2_t[i] = floor(N2_i)

		if t>= t_f:
			break
			
		S.append(S_t)
		N1.append(N1_t)
		N2.append(N2_t)

		

		t += dt
		temps.append(t)
	return S, N1, N2, temps

def main():
	# S_0 = np.array([random.random()*10 for x in range(nb_grad_x)])
	# N1_0 = np.array([random.random()*0 for x in range(nb_grad_x)])
	# N2_0 = np.array([random.random()*10 for x in range(nb_grad_x)])
	# R_S = 1
	# D_N = 5/10000
	# S, N1, N2, temps = resolve_random(S_0, N1_0, N2_0, nu, c, d, D_S, D_N, R_S, R_N1)
	# N1_i, N2_i = integrate(N1), integrate(N2)
	# N1_m, N2_m = moyenne(N1_i[200*int(1/dt):]), moyenne(N2_i[200*int(1/dt):])
	# print("N1_m = {}, N2_m = {}".format(N1_m, N2_m))

	list_R_S = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
	liste_D_N = [25/10000, 10/10000,  5/10000, 2.5/10000]


	l_N1_top = []
	l_N2_top = []

	for d in liste_D_N:
		D_N = d
		liste_int_1 = []
		liste_int_2 = []
		for r in list_R_S:
			R_S = r
			print("D_N = {}, R_S = {}".format(D_N, R_S))
			S_0 = np.array([random.random()*10 for x in range(nb_grad_x)])
			N1_0 = np.array([random.random()*0 for x in range(nb_grad_x)])
			N2_0 = np.array([random.random()*10 for x in range(nb_grad_x)])

			S, N1, N2, temps = resolve_random(S_0, N1_0, N2_0, nu, c, d, D_S, D_N, R_S, R_N1)
			N1_i, N2_i = integrate(N1), integrate(N2)
			N1_m, N2_m = moyenne(N1_i[200*int(1/dt):]), moyenne(N2_i[200*int(1/dt):])
			print("N1_m = {}, N2_m = {}".format(N1_m, N2_m))


			liste_int_1.append(N1_m)
			liste_int_2.append(N2_m+N1_m)

		l_N1_top = liste_int_1 + l_N1_top
		l_N2_top = liste_int_2 + l_N2_top

	fig = plt.figure()
	ax = fig.add_subplot(111, projection="3d")
	# ax.set_ylim([2.5/10000, 5/10000])

	ax.set_yticks(range(len(liste_D_N)))
	ax.set_yticklabels(liste_D_N[::-1])
	ax.invert_yaxis()

	ax.set_xlabel("Ressource influx")
	ax.set_ylabel("Cell diffusion")
	ax.set_zlabel("Population size")
	x_, y_  = np.meshgrid(list_R_S, range(len(liste_D_N)))
	x, y = x_.ravel(), y_.ravel()

	line_red, line_blue = Patch(facecolor='red', edgecolor="black"), Patch(facecolor='blue', edgecolor="black")
	plt.legend([line_blue, line_red], ["N1", "N2"])
	witdh, depth = 0.15, 0.25
	ax.bar3d(x, y, [0]*len(l_N2_top), witdh, depth, l_N1_top, color="blue", shade=True, edgecolor="black")
	ax.bar3d(x, y, l_N1_top, witdh, depth, l_N2_top, color="red", shade=True, edgecolor="black")

	plt.show()



if __name__ == '__main__':
	main()
