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
	def __init__(self, Mtype, Nb, Lmax, Hop, Seed):

		self.Mtype = Mtype
		self.Nb = Nb
		self.Lmax = Lmax
		self.Hop = Hop
		self.Seed = Seed

		self.Lim = 128
		self.Eta = 0.15
		self.Beta = 128
		self.Nquad = 128

		self.model_path = f'./machine_learning/Nb{self.Nb}/{self.Mtype}_noabs'
		self.Create_folder(self.model_path)
		self.model_name = f'{self.Mtype}_Nb{self.Nb}_Hop{self.Hop}_Lmax{self.Lmax}_S{self.Seed}_noabs'	

		self.Nmax, self.Diwn_t, self.Dw_t, self.Bath_s, self.Bath_t, self.Chi_t, self.iwn, self.w, self.Gcff = self.Load_data() #, self.Dcoeff, self.Giwn
		self.Bath_p = self.Model_predict()
		#self.Save_predict()
		#self.Save_optimize()
		#self.Bath_o, self.Chi_o = self.Load_opt()

	def Create_folder(self, folder_path):
		if not os.path.exists(folder_path):
			os.makedirs(folder_path)
			print(f"Folder created: {folder_path}")
		else:
			print(f"Folder already exists: {folder_path}")

	def V_scaling(self, bath):
		bath_s = bath.copy()
		bath_s[:self.Nb] /= np.sqrt( sum( bath_s[:self.Nb]**2 ) )
		#bath_s[:self.Nb] =  bath_s[:self.Nb]
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
			#Giwn = []; Dcoeff = [];
			for n,H in enumerate(self.HGroup):
				group = file[H]	
				diwn = group['Diwn']
				dw = group['Dw']
				gcff = group['Gcoeff']
				bath_set = group[f'bath{self.Nb}/P0']
	
				Bath_t.append(bath_set[1])
				Bath_s.append( self.V_scaling(bath_set[1]) )
	
				Chi_t.append(bath_set.attrs['chi'][0])
					
				Diwn.append( diwn[0] + 1j*diwn[1] );	
		#		Giwn.append( Data['Giwn'][0] + 1j*Data['Giwn'][1] )
					
				Dw.append( dw[0] + 1j*dw[1] )	
					
				Gcff.append( gcff[0][:2*self.Lmax] )
		#		Dcoeff.append( Data['Dcoeff'][0][:2*self.Lmax] )
			if(self.Lmax != 0):
				Gcff = normalize(Gcff, norm="l1")	#normailzation 다른 전처리도 알아보자

		return len(Dw[0]), np.stack(Diwn), np.stack(Dw), np.array(Bath_s), np.array(Bath_t), np.array(Chi_t), 1j*omgn, omg+1j*self.Eta, np.array(Gcff)#, np.array(Dcoeff), np.stack(Giwn)

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
		bath_p = []
		for i,j in enumerate(self.test_list):
			bath_p.append( self.V_Rescaling( bath[i], j) )
	
		return np.array(bath_p)

	def Calculate_hyb(self, bath, omega):
		
		return [sum(bath[l]**2 / (omega[j] - bath[self.Nb + l]) for l in range(self.Nb)) for j in range(self.Nmax)]
	
	def Calculate_Chi(self, bath, index):
		
		Diwn_t = self.Diwn_t[index][:self.Lim]	
		
		Diwn = self.Calculate_hyb(bath, self.iwn)[:self.Lim]
		Diwn = np.array(Diwn)
	
		Chi = np.sum( np.abs(Diwn - Diwn_t)**2 ) / self.Lim
	
		return Chi
		
	def V_Rescaling(self, bath, index):
		dw_t = self.Dw_t[index]
		bath_f = bath.copy()
		bath_f[:self.Nb] /= np.sqrt( sum(bath_f[:self.Nb]**2) )
		bath_f[:self.Nb] *= np.sqrt( -integrate.simpson( np.imag(dw_t), x=np.real(self.w) ) /np.pi )
			
		return bath_f
			
	def Plot_Dw(self, index):
		fig, ax = plt.subplots(2,1, figsize=(12.8,9.6))
		ax[0] = plt.subplot(2,1,1); ax[1] = plt.subplot(2,1,2);
		ax[0].set_xlim(-1.0,.5)
		
		bath_t = self.Bath_t[index]
		dw_t = self.Calculate_hyb(bath_t, self.w)
		ax[0].plot(np.real(self.w), np.imag(dw_t), linewidth=10, color='black', label='true')
		ax[1].scatter(bath_t[self.Nb:],bath_t[:self.Nb], s=300, color='black', label='true')
		
		index_rearange = self.test_list.index(index)

		bath_p = self.Bath_p[index_rearange]
		dw_p = self.Calculate_hyb(bath_p, self.w)
		ax[0].plot(np.real(self.w), np.imag(dw_p), linestyle='--', linewidth=10, label=f'{self.Mtype}$^p$$_{{L={self.Lmax}}}$', alpha=0.5)
		ax[1].scatter(bath_p[self.Nb:], bath_p[:self.Nb], s=300, label=f'${self.Mtype}^p$$_{{L={self.Lmax}}}$', alpha=0.5)
		
	#	bath_o = self.Bath_o[index_rearange]
	#	dw_o = self.Calculate_hyb(bath_o, self.w)
	#	ax[0].plot(np.real(self.w), np.imag(dw_o), label=f'{self.Mtype}$^o$$_{{L={self.Lmax}}}$', alpha=0.5)
	#	ax[1].scatter(bath_o[self.Nb:], bath_o[:self.Nb], label=f'${self.Mtype}^o$$_{{L={self.Lmax}}}$', alpha=0.5)

		v_r2, e_r2, t_r2 = self.R2_score(bath_t, bath_p)
		ax[0].set_title(f'{self.HGroup[index]}', fontsize=30)
		ax[1].set_title(f'R2(V):{v_r2:.2f} R2($\\epsilon$):{e_r2:.2f}, $R2_{{tot}}$:{t_r2:.2f}', fontsize=30)
		ax[0].set_xlabel(f'$\\omega$', fontsize=40)
		ax[0].set_ylabel(f'Im$[\\triangle(\\omega+i\\eta)]$', fontsize=40)
		ax[1].set_xlabel(f'$\\epsilon$', fontsize=40)
		ax[1].set_ylabel(f'$V$', fontsize=40)
		plt.legend(loc='upper right', fontsize=40)
		plt.tight_layout()
		plt.show()

	def Save_predict(self):
		chi_p = []
		for i,j in enumerate(self.test_list):
			chi_p.append( self.Calculate_Chi( self.Bath_p[i], j ) )	
		np.array(chi_p)
		
		save_path = f'{self.model_path}/{self.model_name}.h5'
		if not os.path.exists(save_path):
			with h5py.File(save_path, 'w') as file:
				with h5py.File(f'database/BATH{self.Nb}.h5', 'r') as temp:
					dataset = file.create_dataset('omega', data = temp['omega'])
					for i, group in enumerate(self.HGroup[self.test_list]):
						Group = file.create_group(f'{group}')
						reshaped_test_p = np.reshape(self.Bath_p[i], (1,2*self.Nb))
						dataset = Group.create_dataset('pred_bath', data = reshaped_test_p)
						dataset.attrs['pred_chi'] = chi_p[i]
						Diwn = Group.create_dataset('Diwn', data = temp[f'{group}']['Diwn'])	
		else: 
			print("already exist opt file")
			return;
	
	def Save_optimize(self):
		start_time = time.time()
		print("pred -> C.G. -> optimzation process")
		os.system('make pred_opt')
		os.system(f'./pred_opt {self.Mtype} {self.Nb} {self.Hop} {self.Lmax} {self.Seed}')
		print(f'pred -> optimization done, run_time: {time.time() - start_time}(s)')
	
	def Load_opt(self):
		load_path = f'{self.model_path}/{self.model_name}.h5'

		opt_bath = []; opt_chi = [];	
		with h5py.File(load_path,'r') as file:
			for i, group in enumerate(self.HGroup[self.test_list]):
				true_index = np.where(self.HGroup == group)[0][0]
		
				Group = file[f'{group}']
				pred_bath = Group['pred_bath']
				pred_chi = pred_bath.attrs['pred_chi']
				opt_bath.append( Group['opt_bath'][0] )
				opt_chi.append( Group['opt_bath'].attrs['opt_chi'][0])

		return np.array(opt_bath), np.array(opt_chi)
	
	def R2_score(self, bath_t, bath_p):
		V_r2 = r2_score(bath_t[:self.Nb], bath_p[:self.Nb])
		e_r2 = r2_score(bath_t[self.Nb:], bath_p[self.Nb:])
		tot_r2 = r2_score(bath_t, bath_p)
		return V_r2, e_r2, tot_r2

#		V_r2 = []; e_r2 = []; tot_r2 = [];
#		for i in range(len(self.test_list)):
#			V_r2.append(r2_score(bath_t[i][:self.Nb], bath_p[i][:self.Nb]))
#			e_r2.append(r2_score(bath_t[i][self.Nb:], bath_p[i][self.Nb:]))
#			tot_r2.append(r2_score(bath_t[i], bath_p[i]))
		#	print(f'H: {self.HGroup[self.test_list][i]:12} V:{V_r2[i]:6.4f} ' \
		#			f'e:{e_r2[i]: 6.4f} ' \
		#			f'tot:{tot_r2[i]: 6.4f}')
#		np.array(V_r2); np.array(e_r2); np.array(tot_r2);
#		print(f'mean R2 score> V:{np.mean(V_r2):6.4f} ' \
#				f'e:{np.mean(e_r2): 6.4f} ' \
#				f'tot:{np.mean(tot_r2): 6.4f}')

parser = argparse.ArgumentParser()
parser.add_argument('--model', '-m', type=str, dest='model')
parser.add_argument('--Nb', '-Nb', type=int, dest='Nb')
parser.add_argument('--seed', '-s', type=int, dest='seed')
parser.add_argument('--hop', '-n', type=int, dest='hop') #학습시킬 nearest 최대치를 결정 (5면 1~4까지 모두 학습)
#parser.add_argument('--Lmax', '-L', type=int, dest='Lmax')
args = parser.parse_args();

model = Machine_Learning(Mtype=args.model, Nb=args.Nb, Lmax=5, Hop=args.hop, Seed=args.seed)
model.Plot_Dw(model.test_list[0])
#model.R2_score(model.Bath_t[model.test_list[0]], model.Bath_p[0], model.test_list[0])
