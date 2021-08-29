import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from equation1_odeint import *
"""
Résolution de l'équation 1 (sans terme de diffusion) avec la méthode d'euler (explicite)
1 : respiration
2 : fermentation
"""

def resolve1_euler(initial_values, t_f,  nu=15, d=10, c=1, cste_fct=(32, 30, 2), dt = 1/1000):
	"""initial_values : valeurs des solutions au temps t=0 : (S_0, N1_0, N_20)
	t_f : temps d'arrêt de la simulation
	d,c, nu : paramètres du système d'équa diff
	cste_fct : constantes qui définissent les fonctions (y_1, y_21, y_22)
	dt : pas de temps à utiliser
	"""

	S_0, N1_0, N2_0 = initial_values
	d, c, nu, y_1, y_21, y_22 = parameters

	def maj(S):
		J1_S = S/(1+S)
		J2_S = 100*S/(100+S)
		J1_ATP = y_1*J1_S
		J2_ATP = y_21*J1_S + y_22*J2_S

		return J1_S, J2_S, J1_ATP, J2_ATP


	S = [S_0]
	N1 = [N1_0]
	N2 = [N2_0]
	temps = [0]

	S_t, N1_t, N2_t = S_0, N1_0, N2_0

	t = 0
	while True:
		# print("Temps = ", t)
		J1_S, J2_S, J1_ATP, J2_ATP = maj(S_t)
		S_t = S_t + dt*(nu - N1_t*J1_S - N2_t*J2_S), 
		N1_t =  N1_t*(1 + dt*(c*J1_ATP - d))
		N2_t =  N2_t*(1 + dt*(c*J2_ATP - d))
		if t>t_f:
			break
			
		S.append(S_t)
		N1.append(N1_t)
		N2.append(N2_t)

		t += dt
		temps.append(t)
	return S, N1, N2, temps


def main():
	S_e, N1_e, N2_e, temps_e = resolve_euler(S_0, N1_0, N2_0)
	S_o, N1_o, N2_o, temps = resolve_odeint(S_0, N1_0, N2_0)

	if N1_0 == 0:
		N2_seul = True
	elif N2_0 ==0:
		N1_seul = True

	print(len(temps_e))
	# for  i, t in enumerate(temps_e):
	# 	print(t, " : ", temps[i])

	#Affichage
	fig, ax = plt.subplots()
	ax.set_ylabel("Population size")
	ax.set_xlabel("t")
	ax.tick_params(axis="y")

	ax_S = ax.twinx()
	ax_S.set_ylabel("Ressource size")
	ax_S.tick_params(axis="y", labelcolor="red")
	ax_S.set_ylim([0,1])

	if N1_seul:
		ax_S.plot(temps, S_e[:-2], label='S', color="red")
		ax_S.plot(temps, S_o, '--', color='red')

		ax.plot(temps, N1_e[:-2], label='N1', color = "blue")
		ax.plot(temps, N1_o, '--', color="blue")

	elif N2_seul:
		ax_S.plot(temps, S_e[:-2], label='S', color="red")
		ax_S.plot(temps, S_o, '--', color='red')

		ax.plot(temps, N2_e[:-2], label='N2', color = "orange")
		ax.plot(temps, N2_o, '--', color="orange")

	else:
		ax_S.plot(temps, S_e[:-2], label='S', color="red")
		ax_S.plot(temps, S_o, '--', color='red')

		ax.plot(temps, N1_e[:-2], label='N1', color = "blue")
		ax.plot(temps, N1_o, '--', color="blue")

		ax.plot(temps, N2_e[:-2], label='N2', color = "orange")
		ax.plot(temps, N2_o, '--', color="orange")

	ax.legend()
	ax_S.legend(loc = "center right")
	# first_leg = ax.legend(loc="upper right")
	# plt.gca().add_artist(first_leg)

	# if afficher_cste:
	# 	extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
	# 	plt.legend([extra], [plot_constantes()], loc="upper center")

	plt.show()


if __name__ == '__main__':
	main()

