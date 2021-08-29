import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate
from matplotlib.patches import Rectangle
""" 
Résolution de l'équation 1 (sans terme de diffusion) avec odeint
1 : respiration
2 : fermentation
"""



def resolve1_odeint(initial_values, t_f,  nu=15, d=10, c=1, cste_fct=(32, 30, 2), dt = 1/1000):
	"""initial_values : valeurs des solutions au temps t=0 : (S_0, N1_0, N_20)
	t_f : temps d'arrêt de la simulation
	parameters : paramètres du système d'équa diff : (d, c, nu, y_1, y_21, y_22)
	dt : pas de temps à utiliser
	"""
	S_0, N1_0, N2_0 = initial_values
	y_1, y_21, y_22 = cste_fct

	nb_grad = int(t_f/dt)

	def maj(S):
		J1_S = S/(1+S)
		J2_S = 100*S/(100+S)
		J1_ATP = y_1*J1_S
		J2_ATP = y_21*J1_S + y_22*J2_S
		
		return J1_S, J2_S, J1_ATP, J2_ATP

	def f(y, t):
		S, N1, N2 = y
		J1_S, J2_S, J1_ATP, J2_ATP = maj(S)

		return [nu - J1_S*N1-J2_S*N2, N1*(c*J1_ATP-d), N2*(c*J2_ATP-d)]

	temps = np.linspace(0, t_f, nb_grad)
	sol = scipy.integrate.odeint(f, [S_0, N1_0, N2_0], temps)

	S =  [x[0] for x in sol]
	N1 = [x[1] for x in sol]
	N2 = [x[2] for x in sol]
	return S, N1, N2, temps

if __name__ == '__main__':
	if N1_0 == 0:
		N2_seul = True
	elif N2_0 ==0:
		N1_seul = True

	S, N1, N2, temps = resolve_odeint(S_0, N1_0, N2_0)

	print(S[-10:], N1[-10:], N2[-10:])
	#Affichage
	fig, ax = plt.subplots()
	ax.set_ylabel("Population size")
	ax.set_xlabel("t")
	ax.tick_params(axis="y")

	ax_S = ax.twinx()
	ax_S.set_ylabel("Ressource size")
	ax_S.tick_params(axis="y", labelcolor="green")
	ax_S.set_ylim([0,1])


	if N1_seul:
		ax_S.plot(temps, S, label='S', color="green")
		ax.plot(temps, N1, label='N1', color = "blue")

	elif N2_seul:
		ax_S.plot(temps, S, label='S', color="green")
		ax.plot(temps, N2, label='N2', color = "red")

	else:
		ax_S.plot(temps, S, label='S', color="green")
		ax.plot(temps, N1, label='N1', color = "blue")
		ax.plot(temps, N2, label='N2', color = "red")

	ax.legend()
	ax_S.legend(loc = "center right")

	plt.show()