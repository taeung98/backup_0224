# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('QtAgg')
import math
import h5py
import sys
import os
import re
import time
import numpy as np
import argparse
from joblib import dump, load 
from scipy import integrate
from scipy.special import legendre
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
#from sklearn.model_selection import KFold
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, normalize

class Machine_Learning:
	def __init__(self, Mtype, Nb, Lmax, Hop, Seed, Save=None):
		self.Mtype = Mtype
		self.Nb = Nb
		self.Lmax = Lmax
		self.Hop = Hop
		self.Seed = Seed
	
		self.Lim = 128
		self.Eta = 0.15
		self.Beta = 128
		self.Nquad = 128
	
		self.model_path = f'./machine_learning/Nb{self.Nb}/{self.Mtype}'
		self.Create_folder(self.model_path)
		self.model_name = f'Hop{self.Hop}_Lmax{self.Lmax}_S{self.Seed}'	
		self.Lcoeff_merge()
		stime = time.time()
		self.Nmax, self.Diwn_t, self.Dw_t, self.Bath_s, self.Bath_t, self.Chi_t, self.iwn, self.w, self.Gcff, self.Giwn = self.Load_data()#, self.Dcoeff, self.Giwn
		print(f'Loding data: {time.time() - stime}');

		if Save is not None:
			self.Save_predict()
			self.Save_optimize()
		self.Bath_p, self.Chi_p, self.Bath_o, self.Chi_o = self.Load_pred_opt()

	def Lcoeff_merge(self):
		if not os.path.exists(f'database/BATH{self.Nb}.h5'):
			os.system('make Lcoeff');
			os.system(f'job/job.sh {self.Nb} L')
			os.system(f'py h5/merge.py -Nb {self.Nb}')
		else:
			print(f'BATH{self.Nb}.h5 exist')

	def Create_folder(self, folder_path):
		if not os.path.exists(folder_path):
			os.makedirs(folder_path)
			print(f"Folder created: {folder_path}")
		else:
			print(f"Folder already exists: {folder_path}")

	def V_scaling(self, bath):
		bath_s = bath.copy()
		bath_s[:self.Nb] /= np.sqrt( sum( bath_s[:self.Nb]**2 ) )
		bath_s[:self.Nb] = np.abs( bath_s[:self.Nb] )
		return bath_s

	def Load_data(self):
		with h5py.File(f'database/BATH{self.Nb}.h5', 'r') as file:
			self.HGroup = np.array(list(file.keys()))
			if( 0 < self.Hop < 5 ):
				hop = np.vectorize(lambda x: x.startswith(f'N{self.Hop}'))(self.HGroup)
				self.HGroup = self.HGroup[hop]
			elif( self.Hop > 4 or self.Hop < 0 ):
				print("args.hop erorr: # of nearest neighbor is wrong")
				sys.exit(1);
	
			omg, omgn = file['omega']
			if 'omega' in self.HGroup:
				self.HGroup = np.delete(self.HGroup, np.where( self.HGroup == 'omega' ))
	
			self.train_list, self.test_list = train_test_split( range(len(self.HGroup)), test_size=0.2, random_state=self.Seed )
	
			Diwn = []; Dw = []; Bath_s = []; Chi_t = []; Gcff = []; Bath_t = [];
			Giwn = []; #Dcoeff = [];
			for n,H in enumerate(self.HGroup):
				group = file[H]	
				diwn = group['Diwn']
				dw = group['Dw']
				gcff = group['Gcoeff']
				bath_set = group[f'bath{self.Nb}/P0']
				giwn = group['Giwn']
	
				Bath_t.append(bath_set[1])
				Bath_s.append( self.V_scaling(bath_set[1]) )
	
				Chi_t.append(bath_set.attrs['chi'][0])
				Diwn.append( diwn[0] + 1j*diwn[1] );	
				Giwn.append( giwn[0] + 1j*giwn[1] )
				Dw.append( dw[0] + 1j*dw[1] )	
		#		Dcoeff.append( Data['Dcoeff'][0][:2*self.Lmax] )
					
				if( self.Lmax != 0 ):
					Gcff.append( gcff[0][:2*self.Lmax] )
				else:
					#Gcff.append(diwn[0][:self.Lim]+diwn[1][:self.Lim]);
					Gcff.append(diwn[0]+diwn[1]);
			
			Gcff = normalize(Gcff, norm="l1")	#normailzation
			# Gcff can be thought of as a feature
		return len(Dw[0]), np.stack(Diwn), np.stack(Dw), np.array(Bath_s), np.array(Bath_t), np.array(Chi_t), 1j*omgn, omg+1j*self.Eta, np.array(Gcff), np.stack(Giwn)# np.array(Dcoeff), np.stack(Giwn)

	def Linear_Regression(self):
		model_path = f'{self.model_path}/{self.model_name}.joblib'
		if os.path.exists(model_path):
			model = load(model_path)
		else:
			X_train = self.Gcff[self.train_list]
			y_train = self.Bath_s[self.train_list]
			model = LinearRegression()
			model.fit(X_train, y_train)
			dump(model, model_path)  # model save to .joblib
		return model

	def Model_predict(self):
		if(self.Mtype == 'RF'):
			model = self.Random_Forest()
		elif(self.Mtype == 'LR'):
			model = self.Linear_Regression()
		else:
			print(f'{self.Mtype} model is not exist');
			sys.exit(1);
		
		bath = model.predict(self.Gcff[self.test_list])
		bath_p = []; chi_p = []
		for i,j in enumerate(self.test_list):
			bath_p.append( self.V_Rescaling( bath[i],j) )
			chi_p.append( self.Calculate_Chi( bath_p[i], j) )
	
		return np.array(bath_p), np.array(chi_p)

	def Calculate_hyb(self, bath, omega, ranges=None):
		if ranges is  None:
			ranges = self.Nmax
		
		return np.array([sum(bath[l]**2 / (omega[j] - bath[self.Nb + l]) for l in range(self.Nb)) for j in range(ranges)])

	def Calculate_Chi(self, bath, index):
	
		Chi = np.sum( np.abs( self.Calculate_hyb(bath, self.iwn, self.Lim) - self.Diwn_t[index][:self.Lim] )**2 ) / self.Lim
	
		return Chi

	def V_Rescaling(self, bath, index):
		bath_f = bath.copy()
		bath_f[:self.Nb] *= np.sqrt( -integrate.simpson( np.imag( self.Dw_t[index] ), x=np.real(self.w) ) /np.pi ) / np.sqrt( sum(bath_f[:self.Nb]**2) ) 
			
		return bath_f

	def Plot_Dw(self, index):
#		v_r2, e_r2, t_r2 = self.R2_score(bath_t, bath_p)
		fig, ax = plt.subplots(2,2, figsize=(38.2,20))
		ax[0][0] = plt.subplot(2,2,1); ax[1][0] = plt.subplot(2,2,3);
		ax[0][0].set_xlim(-1.,.5)
		
		colors = ['orange', 'g', 'r']
		alphas = [ 0.5, 0, 0 ]
		fontsizes = [ 30, 35, 40, 45 ]
		line = [ 5, 10, 15, 20, 25, 30, 35 ]

		bath_t = self.Bath_t[index]
		bath_t[:self.Nb] = np.abs(bath_t[:self.Nb])
		dw_t = self.Calculate_hyb(bath_t, self.w)
		ax[0][0].plot(np.real(self.w), np.imag(dw_t), linewidth=line[5], color=colors[0], alpha=alphas[0])
		ax[1][0].scatter(bath_t[self.Nb:], np.abs(bath_t[:self.Nb]), s=2000, marker='x', linewidths=line[0], color=colors[0], alpha=alphas[0])
		re_index = self.test_list.index(index)
	
		bath_p = self.Bath_p[re_index]
		dw_p = self.Calculate_hyb(bath_p, self.w)
		ax[0][0].plot(np.real(self.w), np.imag(dw_p), color=colors[1], linestyle='--', linewidth=line[2])
		ax[1][0].scatter(bath_p[self.Nb:], bath_p[:self.Nb], marker='+', linewidths=line[0], color=colors[1], s=2000)

		bath_o = self.Bath_o[re_index]
		dw_o = self.Calculate_hyb(bath_o, self.w)
		ax[0][0].plot(np.real(self.w), np.imag(dw_o), color=colors[2], linewidth=line[0])
		ax[1][0].scatter(bath_o[self.Nb:], bath_o[:self.Nb], marker='1', linewidths=line[0], color=colors[2], s=2000)

		ax[0][0].set_xlabel(f'$\\omega$', fontsize=fontsizes[2])
		ax[0][0].set_ylabel(f'Im$[\\triangle(\\omega+i\\eta)]$', fontsize=fontsizes[2])
		ax[1][0].set_xlabel(f'$\\epsilon$', fontsize=fontsizes[2])
		ax[1][0].set_ylabel(f'$V$', fontsize=fontsizes[2])

		#imaginary domain
		ax[0][1] = plt.subplot(2,2,2); ax[1][1] = plt.subplot(2,2,4);
		x = np.imag(self.iwn[:self.Lim])

		dw_t = self.Calculate_hyb(bath_t, self.iwn, self.Lim)
		ax[0][1].plot(x, np.real(dw_t), linewidth=line[5], color=colors[0], label='true', alpha=alphas[0])
		ax[1][1].plot(x, np.imag(dw_t), linewidth=line[5], color=colors[0], alpha=alphas[0])

		dw_p = self.Calculate_hyb(bath_p, self.iwn, self.Lim)
		ax[0][1].plot(x, np.real(dw_p), linestyle='--', color=colors[1], linewidth=line[2], label=f'{self.Mtype}$^p$$_{{L={self.Lmax}}}$')#, alpha=alphas[0])
		ax[1][1].plot(x, np.imag(dw_p), linestyle='--', color=colors[1], linewidth=line[2])#, alpha=alphas[0])

		dw_o = self.Calculate_hyb(bath_o, self.iwn, self.Lim)
		ax[0][1].plot(x, np.real(dw_o), color=colors[2], linewidth=line[0], label=f'{self.Mtype}$^o$$_{{L={self.Lmax}}}$')# alpha=0.7)
		ax[1][1].plot(x, np.imag(dw_o), color=colors[2], linewidth=line[0])#, alpha=0.7)

		ax[0][1].set_ylabel(f'Re$[\\triangle(i\\omega_{{n}})]$', fontsize=fontsizes[2])
		ax[1][1].set_xlabel(f'$\\omega_{{n}}$', fontsize=fontsizes[2])
		ax[1][1].set_ylabel(f'Im$[\\triangle(\\omega_{{n}})]$', fontsize=fontsizes[2])

		fig.text(0.9, 0.205, f'$\\chi^{2}_{{t}}$:  {self.Chi_t[index]:.2e}', fontsize=fontsizes[1])
		fig.text(0.9, 0.16, f'$\\chi^{2}_{{p}}$:  {self.Chi_p[re_index]:.2e}', fontsize=fontsizes[1])
		fig.text(0.9, 0.115, f'$\\chi^{2}_{{o}}$:  {self.Chi_o[re_index]:.2e}', fontsize=fontsizes[1])

		ax[0][1].legend(loc='upper right', fontsize=fontsizes[2])

		plt.tight_layout()
#		plt.savefig(f'{self.HGroup[index]}(hyb{self.Nb}_L{self.Lmax})', dpi=300, transparent=True)	
		plt.show()

	def R2_score(self, bath_t, bath_p):
		V_r2 = r2_score(bath_t[:self.Nb], bath_p[:self.Nb])
		e_r2 = r2_score(bath_t[self.Nb:], bath_p[self.Nb:])
		tot_r2 = r2_score(bath_t, bath_p)
		return V_r2, e_r2, tot_r2

	def Save_predict(self):
		path = f'{self.model_path}/{self.model_name}.h5'
		if not os.path.exists(path):
			bath_p, chi_p = self.Model_predict()
			with h5py.File(path, 'w') as file:
				with h5py.File(f'database/BATH{self.Nb}.h5', 'r') as temp:
					dataset = file.create_dataset('omega', data = temp['omega'])
					for i, group in enumerate(self.HGroup[self.test_list]):
						Group = file.create_group(f'{group}')
						reshaped_test_p = np.reshape(bath_p[i], (1,2*self.Nb))
						dataset = Group.create_dataset('pred_bath', data = reshaped_test_p)
						dataset.attrs['pred_chi'] = chi_p[i]
						Diwn = Group.create_dataset('Diwn', data = temp[f'{group}']['Diwn'])	

	def Save_optimize(self):
		start_time = time.time()
		print("pred -> C.G. -> optimzation process")
		os.system('make pred_opt')
		os.system(f'./pred_opt {self.Mtype} {self.Nb} {self.Hop} {self.Lmax} {self.Seed}')
		print(f'pred -> optimization done, run_time: {time.time() - start_time}(s)')

	def Load_pred_opt(self):
		load_path = f'{self.model_path}/{self.model_name}.h5'
		bath_p = []; chi_p = []; opt_bath = []; opt_chi = [];	
		with h5py.File(load_path,'r') as file:
			for i, group in enumerate(self.HGroup[self.test_list]):
				true_index = np.where(self.HGroup == group)[0][0]
				Group = file[f'{group}']
				pred_bath = Group['pred_bath']
				bath_p.append( pred_bath[0] )
				chi_p.append(pred_bath.attrs['pred_chi'])
				opt_bath.append( Group['opt_bath'][0] )
				opt_chi.append( Group['opt_bath'].attrs['opt_chi'][0])
		return np.array(bath_p), np.array(chi_p), np.array(opt_bath), np.array(opt_chi)

parser = argparse.ArgumentParser();
parser.add_argument('--model', '-m', type=str, dest='model')
parser.add_argument('--Nb', '-Nb', type=int, dest='Nb')
parser.add_argument('--seed', '-s', type=int, dest='seed')
parser.add_argument('--hop', '-n', type=int, dest='hop') #학습시킬 nearest 최대치를 결정 (5면 1~4까지 모두 학습)
parser.add_argument('--Lmax', '-L', type=int, dest='Lmax')
args = parser.parse_args();

try:

	model = Machine_Learning(Mtype=args.model, Nb=args.Nb, Lmax=args.Lmax, Hop=args.hop, Seed=args.seed)
	model.Plot_Dw(model.test_list[0])
#	model.Plot_Diwn(model.test_list[0])

except KeyboardInterrupt:
	print("")
print("END")
