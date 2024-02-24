# ml.py class Machine_Learing 내장 함수 빼놈	 

def FT(self, indx=0):
	ti, wi = np.polynomial.legendre.leggauss(self.nquad)	
	taus = 0.5*(ti + 1) * self.Beta
	Gtau = np.array(
		[np.sum([ np.exp(-w*t) * g for w, g in zip(self.iwn, self.Giwn[indx] - 1/self.iwn)]) for t in taus])
	return taus, Gtau/self.Beta - 0.5, wi

def LT(self, index=0):
	self.tau, self.Gtau, self.tw = self.FT(indx=index)
	Gcoeffs = []; 
	for p in range(self.Lmax):
		add = 0
		for t, g, w in zip(self.tau, self.Gtau.real, self.tw):
			add += g * w * legendre(p)(2*t/self.Beta-1)
		Gcoeffs += [ add *(2*p+1)/2 ]
	result = [ val for pair in zip(np.real(Gcoeffs), np.imag(Gcoeffs)) for val in pair ]
	print(index)	
	return result

def ILT(self):
	return np.sum( [ self.Gcoeff[p] * legendre(p)(2*self.tau/self.beta-1) for p in range(0, self.Lmax, 2) ], axis=0)
	
def Show_LT(self):
	plt.scatter(self.tau, self.ILT(), label='Legendre')
	plt.plot(self.tau, self.Gtau.real, label='Fourier', color='r')
	plt.show()
