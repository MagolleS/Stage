import numpy as np 
from fonctions import *
import random
import pickle
from math import floor, ceil

""" 
Résolution de l'équation 2  en 2D, avec la méthode d'euler pour la variable temporelle et odeint pour la variable spatiale.
1 : respiration
2 : fermentation
Les matrices M_inv_S et M_inv_N sont très lourdes, pour éviter de les calculer à chaque fois, 
je les stocke avec pickle dans n fichier qui fait environ 1.5 Go.
"""
x_0, x_f = 0,1
y_0, y_f = 0,1
def resolve2_euler_2d(initial_values, t_f, t_N1_0, t_N1_f, nu=15, d=10, c=1, D_N=5/10000, D_S=0.25/10000, R_S=0, R_N1=0, proba_ressource=0, proba_N1=0, arrondir_cellules=False, cste_fct=(32, 30, 2), dt = 1/1000, dx = 1/100, dy=1/100):

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
	nb_grad_y = int((y_f - y_0)/dy)

	n = nb_grad_x*nb_grad_y

	nu_l = np.array([nu]*n)
	d_l = np.array([d]*n)
	c_l = np.array([c]*n)


	S, N1, N2 = [S_0], [N1_0], [N2_0]
	S_t, N1_t, N2_t = S_0, N1_0, N2_0

	def maj(S):
		J1_S = S/(1+S)
		J2_S = 100*S/(100+S)
		J1_ATP = y_1*J1_S
		J2_ATP = y_21*J1_S + y_22*J2_S

		return J1_S, J2_S, J1_ATP, J2_ATP



	# A = laplacien2DPeriodic(nb_grad_x, nb_grad_y)	
	# M_inv_S = np.linalg.inv((np.identity(nb_grad_x*nb_grad_y)+D_S*dt*A))
	# M_inv_N = np.linalg.inv((np.identity(nb_grad_x*nb_grad_y)+D_N*dt*A))


	# with open("matrice_save.obj", 'wb') as f:
	# 	pickle.dump((M_inv_S, M_inv_N), f)

	with open("matrice_save.obj", 'rb') as f:
		M_inv_S, M_inv_N = pickle.load(f)


	temps = [0]
	t = 0
	while True:
		if abs(t-int(t))<10**(-3):
			print(t)
		J1_S, J2_S, J1_ATP, J2_ATP = maj(S_t)

		S_influx = np.array([0]*n)
		for i in range(n):
			if random.random() < proba_ressource:
				S_influx[i] += R_S

		N1_influx = np.array([0]*n)	
		if t_N1_0 <= t <= t_N1_f:
			for i in range(n):
				if random.random() < proba_N1:
					N1_influx[i] += R_N1

		N2_influx = np.array([0]*n)

		S_t = M_inv_S.dot(S_t + dt*(nu_l - N1_t*J1_S - N2_t*J2_S)) + S_influx
		N1_t = M_inv_N.dot(N1_t + dt*N1_t*(c_l*J1_ATP - d_l)) + N1_influx
		N2_t = M_inv_N.dot(N2_t + dt*N2_t*(c_l*J2_ATP - d_l)) + N2_influx

		if arrondir_cellules:
			for i in range(n):
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
	# Affichage
	afficher2D(S, N1, N2, temps)
	# afficherIntegrate(S,N1, N2, temps, xy=True)

if __name__ == '__main__':
	main()
