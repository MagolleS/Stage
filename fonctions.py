import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm
import animatplot as amp


def laplacienPeriodic(n, correction=True):
	#renvoie la matrice du laplacien avec condition périodique au bord en 1 dimension, de taille n*n
	A = 2*np.identity(n)
	for i in range(n):
		if i < n-1:
			A[i][i+1] = -1
		if i>0:
			A[i][i-1] = -1

	A[0][n-1] = -1
	A[n-1][0] = -1
	if correction:
		return A

	return -A

def laplacien2DPeriodic(n, m):
	#renvoie la matrice du laplacien avec condition périodique au bord en 2 dimension, de taille (n*m)x(n*m)

	i_dx2 = -1/dx2
	i_dy2 = -1/dy2
	A = -2*(i_dx2+i_dy2)*np.identity(n*m)
	for k in range(m):
		for j in range(n):
			if j == 0:
				A[n*k+j][n*(k+1)-1] = i_dx2
			else:
				A[n*k+j][n*k+j-1] = i_dx2

			if j==n-1:
				A[n*k+j][n*k] = i_dx2
			else:
				A[n*k+j][n*k+j+1] = i_dx2

			if k==0:
				A[n*k+j][n*(m-1) + j] = i_dy2
			else:
				A[n*k+j][n*(k-1) + j] = i_dy2

			if k == m-1:
				A[n*k+j][j] = i_dy2
			else:
				A[n*k+j][n*(k+1) + j] = i_dy2

	return A

def integrateX(N, dx):
	return dx * sum(N)

def integrateXY(N, dx, dy):
	return dx*dy*sum(N)

def integrate2D(S, dx, dy):
	return [integrateXY(s, dx, dy) for s in S]

def integrate(S, dx):
	return [integrateX(s, dx) for s in S]

def moyenne(l):
	return sum(l)/len(l)


def afficher(S, N1, N2, temps, show_integrate=True):
	fig, ax = plt.subplots()
	ax.set_ylabel("Population size")
	ax.set_xlabel("x")
	ax.tick_params(axis="y")
	lineN1, = ax.plot(espaceX, N1[0], label="N1", color = "blue")
	lineN2, = ax.plot(espaceX, N2[0], label="N2", color = "orange")

	ax_S = ax.twinx()
	ax_S.set_ylabel("Ressource size")
	ax_S.tick_params(axis="y", labelcolor="green")
	lineS, = ax_S.plot(espaceX, S[0], label="S", color="green")

	text = plt.text(0.91,0.77,'t = {}'.format(temps[0]) ,horizontalalignment='center',
		     verticalalignment='center', transform = ax.transAxes)

	def animate(i):
		i*=100
		y_max = max(np.max(N1[i]), np.max(N2[i]))

		ax.set_ylim([0, max(50	,1.05*y_max)])
		ax_S.set_ylim([0, max(1, 1.05*np.max(S[i]))])

		lineS.set_ydata(S[i])
		lineN1.set_ydata(N1[i])
		lineN2.set_ydata(N2[i])
		text.set_text("t = " + str(temps[i])[:6])

	anim = FuncAnimation(fig, animate, interval=dt, frames = len(temps)//100)

	fig.tight_layout()
	ax.legend()
	ax_S.legend()
	plt.draw()
	plt.show()



	if show_integrate:
		afficherIntegrate(S, N1, N2, temps)
		

def afficherIntegrate(S, N1, N2, temps, xy=False):
	if xy:
		S = integrate2D(S, 1/100, 1/100)
		N1 = integrate2D(N1, 1/100, 1/100)
		N2 = integrate2D(N2, 1/100, 1/100)
	else:
		S = integrate(S, 1/100)
		N1 = integrate(N1, 1/100)
		N2 = integrate(N2, 1/100)

	fig, ax = plt.subplots()
	ax.set_ylabel("Population size")
	ax.set_xlabel("t")
	ax.tick_params(axis="y")
	y_max = max(np.max(N1), np.max(N2))
	ax.set_ylim([0, 1.05*y_max])

	ax_S = ax.twinx()
	ax_S.set_ylabel("Ressource size")
	ax_S.tick_params(axis="y", labelcolor="green")
	ax_S.set_ylim([0, 1.05*np.max(S)])

	ax_S.plot(temps, S, label='S', color="green")
	ax.plot(temps, N1, label='N1', color = "blue")
	ax.plot(temps, N2, label='N2', color = "orange")

	ax_S.legend(loc="center right")
	ax.legend(loc="upper right")
	plt.show()


	print("S : ", S[-10:])
	print("N1 : ", N1[-10:])
	print("N2 : ", N2[-10:])


def afficherSurface(S, N1, N2, temps):
	X, T = np.meshgrid(espaceX, temps)
	plt.xlabel('x')
	plt.ylabel('t', rotation=0)

	plt.pcolormesh(X, T, S, shading="auto")
	plt.colorbar()
	plt.title("S")
	plt.show()

	plt.xlabel('x')
	plt.ylabel('t', rotation=0)

	plt.pcolormesh(X, T, N1, shading="auto")
	plt.colorbar()
	plt.title("N1")
	plt.show()


	plt.xlabel('x')
	plt.ylabel('t', rotation=0)
	
	plt.pcolormesh(X, T, N2, shading="auto")
	plt.colorbar()
	plt.title("N2")
	plt.show()




def afficherBarre(S, N1, N2, temps):
	S = integrate(S)
	N1 = integrate(N1)
	N2 = integrate(N2)

	xs = np.array([1, 2, 3])

	fig, ax = plt.subplots()

	ax.set_xticks(xs)

	barreN1 = ax.bar(2, [N1[0]], 0.5, label = "N1", align="center", color="blue")
	barreN2 = ax.bar(2, [N2[0]], 0.5, bottom=[N1[0]], label = "N2", align="center", color="orange")

	ax.bar(1, [0])
	ax.bar(3, [0])

	text = plt.text(0.91,0.77,'t = {}'.format(temps[0]) ,horizontalalignment='center',
		     verticalalignment='center', transform = ax.transAxes)

	def animate(i):
		ax.set_ylim([0, max(50, 1.05*(N1[i] + N2[i
			]))])

		text.set_text("t = " + str(temps[i])[:6])

		barreN1.patches[0].set_height(N1[i])
		barreN2.patches[0].set_y(N1[i])
		barreN2.patches[0].set_height(N2[i])

		# barreN1 = ax.bar(2, [N1[i]], 0.5, label = "N1", align="center", color="blue")
		# barreN2 = ax.bar(2, [N2[i]], 0.5, bottom=[N1[i]], label = "N2", align="center", color="orange")



	anim = FuncAnimation(fig, animate, interval=dt, frames = len(temps))

	plt.legend()
	plt.draw()
	plt.show()
def vector_to_matrix(A, n, m):
	return [A[i*n:(i+1)*n] for i in range(m)]

def afficher2D(S_, N1_, N2_, temps, nb_grad_x, nb_grad_y, show_integrate=True):
	espaceX = np.linspace(0, 1, nb_grad_x)
	espaceY = np.linspace(0, 1, nb_grad_y)
	x, y, t = np.meshgrid(espaceX, espaceY, temps)


	S, N1, N2 = [], [], []
	for i in range(len(temps)):
		S.append(vector_to_matrix(S_[i], nb_grad_x, nb_grad_y))
		N1.append(vector_to_matrix(N1_[i], nb_grad_x, nb_grad_y))
		N2.append(vector_to_matrix(N2_[i], nb_grad_x, nb_grad_y))

	plt.figure()
	plt.title('Ressource')
	y_maxS = max([max([np.max(s_x) for s_x in s_y]) for s_y in S])
	print("y_maxS = ", y_maxS)
	blockS = amp.blocks.Pcolormesh(y[:,:,0], x[:,:,0], S, vmin=0, vmax=y_maxS)
	plt.colorbar(blockS.quad)

	timeline = amp.Timeline(t, units='s', fps=20)
	animS = amp.Animation([blockS], amp.Timeline(t))
	animS.controls()

	plt.figure()
	plt.title("Respirateurs")
	y_maxN1 = max([max([np.max(n1_x) for n1_x in n1_y]) for n1_y in N1])
	print("y_maxN1 = ", y_maxN1)
	blockN1 = amp.blocks.Pcolormesh(y[:,:,0], x[:,:,0], N1, vmin=0, vmax=y_maxN1)
	plt.colorbar(blockN1.quad)

	animN1 = amp.Animation([blockN1], amp.Timeline(t))
	animN1.controls()

	plt.figure()
	plt.title("Fermentateurs")
	y_maxN2 = max([max([np.max(n2_x) for n2_x in n2_y]) for n2_y in N2])
	print("y_maxN2 = ", y_maxN2)
	blockN2 = amp.blocks.Pcolormesh(y[:,:,0], x[:,:,0], N2, vmin=0, vmax=y_maxN2)
	plt.colorbar(blockN2.quad)

	animN2 = amp.Animation([blockN2], amp.Timeline(t))
	animN2.controls()
	plt.show()

	if show_integrate:
		afficherIntegrate(S_,N1_, N2_, temps, xy=True)



def afficherTest(S, N1, N2, temps, show_integrate=True):
	fig, ax = plt.subplots()
	ax.set_xlabel("x")
	ax.set_ylabel("Population size")

	y_max = max(max([np.max(N1i) for N1i in N1]), max([np.max(N2i) for N2i in N2]))
	ax.set_ylim([0, 1.05*y_max])

	axS = ax.twinx()
	axS.set_ylabel("Ressource size")

	y_max = max([np.max(s) for s in S])
	axS.set_ylim([0, 1.05*y_max])


	x, t = np.meshgrid(espaceX, temps)
	lineS = amp.blocks.Line(x,S, ax=axS, color="green", label="S")
	lineN1 = amp.blocks.Line(x,N1, ax=ax, label="N1", color = "blue")
	lineN2 = amp.blocks.Line(x,N2, ax=ax, label="N2", color = "orange")


	timeline = amp.Timeline(t, units='s', fps=20)
	anim = amp.Animation([lineS, lineN1, lineN2], timeline)

	anim.controls()
	ax.legend()
	axS.legend()
	plt.show()

	if show_integrate:
		afficherIntegrate(S, N1, N2, temps)

