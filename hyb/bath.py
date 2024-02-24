import h5py 
import matplotlib.pyplot as plt
from matplotlib import gridspec
import argparse
import numpy as np
from matplotlib.ticker import FuncFormatter
import math

'''	
def extract_variables(self, string):
	var_dict = {}
	pattern = re.compile(r'[-]?\d*\.\d+|\d+')
	sp = re.split(r'[/_]', string)
	for var in sp:
		num = re.findall(pattern, var)
		strings = [_f for _f in re.split(pattern, var) if _f]
		var_dict.update(dict(list(zip(strings, num))))
	return var_dict
'''

class BATH_FITTING:
	def __init__(self, Bnum, Order=False, tend=False):
		self.Nb = Bnum;
		self.O = Order;
		self.snum, self.bath, self.chi, self.wnum, self.omega, self.D, self.G = self.call_data()
		self.tend = tend
		self.ALL_H = [];
	#	self.file_name = 'BATH_N1.h5'

	def format_func(self, value, tick_number):
		return f'{value:.2e}';

	def call_data(self):
		snum = 10;
		bath = [];	chi = [];
		D = [[], []]; G = [[], []];
		file_name = 'test.h5'
		with h5py.File( file_name, 'r') as file:
			self.HGroup = np.array(list(file.keys()))
			if 'omega' in self.HGroup:
				self.HGroup = np.delete(self.HGroup, np.where( self.HGroup == 'omega' ))
			self.H = self.HGroup[0]
			omega = [file[f'/omega']];
			D[0] = [ r + i*1j for r,i in zip([file[f'{self.H}/Dw']][0][0], [file[f'{self.H}/Dw']][0][1])];	
			D[1] = [ r + i*1j for r,i in zip([file[f'{self.H}/Diwn']][0][0], [file[f'{self.H}/Diwn']][0][1])];	
#			G[0] = [ r + i*1j for r,i in zip([file[f'{self.H}/Gw']][0][0], [file[f'{self.H}/Gw']][0][1])];	
#			G[1] = [ r + i*1j for r,i in zip([file[f'{self.H}/Giwn']][0][0], [file[f'{self.H}/Giwn']][0][1])];	

			wnum = omega[0].attrs['Nmax'][0];
			omega = np.array(omega[0]);	
			bath = [file[f'{self.H}/bath{self.Nb}/P0']];
			chi_arr = bath[0].attrs['chi'];
			chi.append(math.floor(chi_arr[0]*1e14)*1e-14);
			for i in range(1,snum):
				add_data = file[f'{self.H}/bath{self.Nb}/P{i}'];
				chi_arr = np.vstack( (chi_arr, add_data.attrs['chi']) );
				bath = np.vstack( (bath, [add_data]) );
				chi.append(math.floor(chi_arr[i][0]*1e14)*1e-14);
				
		return snum, bath, chi, wnum, omega, D, G;

	def call_hyb(self, domain, bath):
		ETA = 0.15;	#brodening rate도 소스코드에서 받아오도록
		hyb = [[], []];
		if domain=='r': 
			for i in range(self.wnum):
				delta1 = 0.
				delta2 = 0.
				for j in range(self.Nb):
					delta1 += bath[0][j]**2/(self.omega[0][i] + 1j*ETA - bath[0][self.Nb + j]) # inital
					delta2 += bath[1][j]**2/(self.omega[0][i] + 1j*ETA - bath[1][self.Nb + j]) # final
				hyb[0].append(delta1)
				hyb[1].append(delta2)
				
		else: 
			for i in range(self.wnum):
				delta1 = 0.
				delta2 = 0.
				for j in range(self.Nb):
					delta1 += bath[0][j]**2/(self.omega[1][i] * 1j - bath[0][self.Nb + j])
					delta2 += bath[1][j]**2/(self.omega[1][i] * 1j - bath[1][self.Nb + j])
				hyb[0].append(delta1)
				hyb[1].append(delta2)
		return hyb
	
	def plot_bath(self):
		fig = plt.figure(figsize=(21,9));	
		gs = gridspec.GridSpec(nrows=1, ncols=2, height_ratios=[6], width_ratios=[12, 1])
		snum, bath, chi, Nb = self.snum, self.bath, self.chi, self.Nb;
		axis = plt.subplot(gs[0]);	
		axis.set_xlim(-.5,.5)
		colors = ['blue', 'red'];
		labels = ['Initial', 'Final'];
		for i in range(snum):
			y_ticks = np.full(Nb, i);
			for j in range(2):
				axis.scatter(bath[i][j][Nb:], y_ticks+.2*j, color=colors[j], marker='o', s=20, label=labels[j] if i==0 else None);
				for k in range(Nb):
					#              ep-V/2 ~ ep+V/2
					axis.plot([bath[i][j][Nb+k]-bath[i][j][k]/2, bath[i][j][Nb+k]+bath[i][j][k]/2], [i+.2*j, i+.2*j], color=colors[j], lw=3);
					axis.scatter([bath[i][j][Nb+k]-bath[i][j][k]/2, bath[i][j][Nb+k]+bath[i][j][k]/2], [i+.2*j, i+.2*j], color=colors[j], marker='|', s=50);
		
		axis.set_yticks([]);
		axis.set_xlabel('$\epsilon$', fontsize=40);
		axis.set_ylabel('Index of initial conditions', fontsize=30);
		axis.legend(bbox_to_anchor=(1., 1), loc='upper left', fontsize=20)
#		chi plot 
#		chi = self.call_data()[2];
		ax = plt.subplot(gs[1], sharey=plt.subplot(gs[0]) );	
		for i in range(snum):
			ax.axvline(chi[i], color='gray', linestyle=':');
			ax.scatter(chi[i], 0.1 + i, color='green', marker='o', s=20);
	
		ax.xaxis.set_major_formatter(FuncFormatter(self.format_func));
		ax.set_xticks([chi[0], chi[-1]])
		plt.tick_params(axis='x', length=2, width=2, labelsize=15)
		ax.set_xlabel('$\chi^2$', fontsize=25);	
		ax.xaxis.set_label_coords(.4, -.033);
#		ax.semilogx(base=10)
		plt.suptitle(f'$N_b={Nb}$', fontsize=30)
		plt.tight_layout()
		plt.savefig(f'{self.H}bath_N{Nb}.png', dpi=500, transparent=True);

	def plot_DOS(self):
		fig = plt.figure()
		plt.plot(self.omega[0][3584:4608], np.imag(self.D[0][3584:4608])*(-1/np.pi))
		hyb = [self.call_hyb('r',self.bath[0]), self.call_hyb('i',self.bath[0])];
		plt.plot(self.omega[0][3584:4608], np.imag(hyb[0][1][3584:4608])*(-1/np.pi))
		plt.ylabel("$\Delta(\omega)$", fontsize = 20)
		plt.xlabel("$\omega$",fontsize = 20)
		plt.title(f"N_{b}={self.Nb}", fontsize = 20);
		plt.xticks([])
		plt.yticks([]); plt.ylim(0,.15)
		plt.savefig(f'hyb{self.Nb}.png', dpi=500, transparent=True);
	
	def subplot_hyb(self, hyb, domain, part):
		if domain=='r':
			if part=='real':
				plt.scatter(self.omega[0], np.real(self.D[0]),  s = 8, label="o")
				plt.scatter(self.omega[0], np.real(hyb[0][1]),  s = 3, label="f")
			else:
				plt.scatter(self.omega[0], np.imag(self.D[0]),  s = 8, label="o")
				plt.scatter(self.omega[0], np.imag(hyb[0][1]),  s = 3, label="f")
		else:
			if part=='real':
				plt.scatter(self.omega[1][0:50], np.real(self.D[1][0:50]),  s = 8, label="o")
				plt.scatter(self.omega[1][0:50], np.real(hyb[1][1][0:50]),  s = 3, label="f")
			else:
				plt.scatter(self.omega[1][0:50], np.imag(self.D[1][0:50]),  s = 8, label="o")
				plt.scatter(self.omega[1][0:50], np.imag(hyb[1][1][0:50]),  s = 3, label="f")

	def plot_hyb(self):
		order = self.O-1;
		hyb = [self.call_hyb('r',self.bath[order]), self.call_hyb('i',self.bath[order])];
		plt.figure()	
		for i,j in enumerate(['real','imag']):
			plt.subplot(2, 2, i+1)
			self.subplot_hyb(hyb, 'r', j);
			plt.subplot(2, 2, i+3)
			self.subplot_hyb(hyb, 'i', j);
		plt.legend(bbox_to_anchor=(1., 1), loc='upper left');

	def overlap_hyb(self):
		plt.figure()
		plt.plot(self.omega[0], np.imag(self.D[0]), color='black', label='Target')
		for i in range(self.snum):
			plt.plot(self.omega[0], np.imag(self.call_hyb('r', self.bath[i])[1]), linestyle='--', label=f'{i}')
#		plt.ylabel("DOS")
		plt.legend();

	def store_top_groups(self, name, obj):
		if isinstance(obj, h5py.Group) and "/" not in name:
			self.ALL_H.append(name)

	def	plot_tend(self):
#		Nb = list(range(3, 21, 2));
		file_name = 'BATH5.h5'
		with h5py.File( file_name, 'r') as file:
			file.visititems(self.store_top_groups)
			for H in self.ALL_H:
				Nb_list = []
				fig = plt.figure();
				group1 = file[H]
				for group_name in group1:
					if group_name.startswith('bath'):
						try:
							Nb = int(group_name[4:])
							Nb_list.append(Nb)
						except ValueError:
							print(f"Invalid number format in group name: {group_name}")
				chi_arr = [];
				print(Nb_list)
				for n in Nb_list: 
					chi = math.floor( file[f'{H}/bath{n}/P0'].attrs['chi'] *1e15 ) * 1e-15
					chi_arr.append(chi)
				plt.scatter( Nb_list, chi_arr );

				plt.semilogy(base=10)
				plt.xticks(Nb_list)
				plt.ylabel('$\chi^2$',fontsize=20 )
				plt.tick_params(axis='both', which='major', labelsize=12)
				plt.ylabel("${\chi}^2$", fontsize = 20)
				plt.xlabel("$N_b$",fontsize = 20)
		#		plt.title("$H_{latt}$", fontsize = 20);
				plt.tight_layout()
				plt.savefig(f'Nb-chi{self.H}.png', dpi=500, transparent=True);

	def total_plot(self):
		try:
			self.plot_bath()
			if self.O: self.plot_hyb()  
			if self.tend: self.plot_tend()
#			self.overlap_hyb()	
#			self.plot_DOS()
			plt.show()
		except KeyboardInterrupt:
			# Handle Ctrl+C, if needed
			print(""); print("Plotting interrupted by user.")
	
# os.path.join: connect to path, os.getcwd(): get current working directory
parser = argparse.ArgumentParser();
#parser.add_argument('--hamiltonian', '-H', type=str, dest='hamiltonian');
parser.add_argument('--bathnum', '-Nb', type=int, nargs='+', dest='bathnum');
parser.add_argument('--order', '-O', type=int, nargs='+', dest='order');
parser.add_argument('--tendency', '-t', type=int, nargs='+', dest='tendency');
args = parser.parse_args();

Nb = args.bathnum[0];	

plot_generator = BATH_FITTING(Nb)#, args.order[0], args.tendency[0]);
plot_generator.total_plot()

