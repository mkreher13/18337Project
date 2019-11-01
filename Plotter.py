#18.337 Project Miriam Kreher
#Class to plot geometry and rays

import numpy as np
import matplotlib.pyplot as plt

class Plotter():

	def __init__(self):
		self

######################################################

	def pincell(self, Pitch, rad):
		HalfPitch = Pitch/2.

		theta = np.linspace(0, 2*np.pi, 100)

		fig, self.ax = plt.subplots()
		self.ax.plot([-HalfPitch,-HalfPitch],[-HalfPitch,HalfPitch], color='blue')
		self.ax.plot([-HalfPitch,HalfPitch],[HalfPitch,HalfPitch], color='blue')
		self.ax.plot([HalfPitch,HalfPitch],[HalfPitch,-HalfPitch], color='blue')
		self.ax.plot([HalfPitch,-HalfPitch],[-HalfPitch,-HalfPitch], color='blue')

		for r in rad:
			x = r*np.cos(theta)
			y = r*np.sin(theta)
			self.ax.plot(x,y)
		self.ax.set_aspect(1)

######################################################

	def rays(self, X1, Y1, X2, Y2):

		self.ax.plot([X1,X2],[Y1,Y2], color='red')


	def end(self):
		plt.savefig("pincell")
		



