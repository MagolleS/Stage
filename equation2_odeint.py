import numpy as np 
import scipy.integrate
from fonctions import *

""" 
Résolution du système d'équa diff avec diffusion en 1D
"""
x_0, x_f = 0,1
def resolve2_odeint(initial_values, t_f, nu=15, d=10, c=1, D_N=5/10000, D_S=0.25/10000, cste_fct=(32, 30, 2), dt = 1/1000, dx = 1/100):

	"""initial_values : valeurs des solutions au temps t=0 : (S_0, N1_0, N_20)
	t_f : temps d'arrêt de la simulation
	d,c, nu, D_N, D_S : paramètres du système d'équa diff
	cste_fct : constantes qui définissent les fonctions (y_1, y_21, y_22)
	dt : pas de temps à utiliser
	"""
	S_0, N1_0, N2_0 = initial_values
	y_1, y_21, y_22 = cste_fct

	nb_grad_t = int(t_f/dt)
	nb_grad_x = int((x_f - x_0)/dx)
	dx2 = dx**2

	def maj(S):
		J1_S = S/(1+S)
		J2_S = 100*S/(100+S)
		J1_ATP = y_1*J1_S
		J2_ATP = y_21*J1_S + y_22*J2_S

		return J1_S, J2_S, J1_ATP, J2_ATP

	nu_l = np.array([nu]*nb_grad_x)
	d_l = np.array([d]*nb_grad_x)
	c_l = np.array([c]*nb_grad_x)

	espace = np.linspace(x_0, x_f, nb_grad_x)
	temps = np.linspace(0, t_f, nb_grad_t)

	A = laplacienPeriodic(nb_grad_x)
	M_S = D_S/dx2 * A
	M_N = D_N/dx2 * A

	def f(X, t):
		S = X[:nb_grad_x]
		N1 = X[nb_grad_x:2*nb_grad_x]
		N2 = X[2*nb_grad_x:]

		J1_S, J2_S, J1_ATP, J2_ATP = maj(S)
		
		S_t = -M_S.dot(S) + nu_l - J1_S*N1 - J2_S*N2
		N1_t = -M_N.dot(N1) + N1*(c_l*J1_ATP - d_l)
		N2_t = -M_N.dot(N2) + N2*(c_l*J2_ATP - d_l)

		return np.concatenate([S_t, N1_t, N2_t])


	sol = scipy.integrate.odeint(f, np.concatenate([S_0, N1_0, N2_0]), temps)

	S = [x[:nb_grad_x] for x in sol]
	N1 =[x[nb_grad_x:2*nb_grad_x] for x in sol]
	N2 = [x[2*nb_grad_x:] for x in sol]

	return S, N1, N2, temps
	#Affichage
	# afficher(S, N1, N2, temps)