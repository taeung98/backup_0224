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
import subprocess

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
	def __init__(self, Mtype, Nb, Lmax, Hop, random_state):#, n_estimate=False, depth=None):
		self.Mtype = Mtype
#		self.plot_already_done = Plot_option
		self.Nb = Nb
		self.lim = 128
		self.eta = 0.15
		self.beta = 128
		self.Lmax = Lmax
		self.nquad = 128
		self.Hop = Hop
		self.seed = random_state  
#		self.n_estimate = n_estimate
#		self.dmax = depth
		self.model_path = f'./machine_learning/Nb{self.Nb}/{self.Mtype}'
		self.create_folder_if_not_exists(self.model_path)
		self.model_name = f'{self.Mtype}_Nb{self.Nb}_Hop{self.Hop}_Lmax{self.Lmax}_S{self.seed}.joblib'	
		
		self.Nmax, self.Diwn, self.Dw, self.bath, self.bath_true, self.chi, self.iwn, self.w, self.Gcoeff = self.Load_data() #, self.Dcoeff, self.Giwn
	
		self.model = self.Select_model()
		self.train_p, self.test_p = self.Predict()
		self.chi_p, self.Diwn_p = self.Chi(self.test_p)
		self.Save_pred()
		self.opt_bath, self.opt_chi = self.Pred_opt()
	
		self.opt_chi, self.opt_Diwn = self.Chi(self.opt_bath)
	
	def Claculate_hyb(self, bath):
		Diwn = []
		Diwn.append([sum(bath[l]**2 / (self.iwn[j] - bath[self.Nb + l]) for l in range(self.Nb)) for j in range(self.lim)])
		Diwn = np.array(Diwn)
	
		return Diwn

	def Re_normalize(self, bath_i):
	
		Dw_r = self.Calculate_hyb(bath_i)
		bath_f = bath_i.copy()
		bath_f[:self.Nb] *= np.sqrt( -integrate.simpson( np.imag(Dw_r[i]), x=self.w ) /np.pi )
	
		return bath_f
	
	def create_folder_if_not_exists(self, folder_path):
		if not os.path.exists(folder_path):
			os.makedirs(folder_path)
			print(f"Folder created: {folder_path}")
		else:
			print(f"Folder already exists: {folder_path}")

	def Select_model(self):
		if(self.Mtype == 'RF'):
			return self.Random_Forest()
		elif(self.Mtype == 'LR'):
			return self.Linear_Regression()
		else:
			print(f'{self.Mtype} model is not exist');
			sys.exit(1);

	def Load_data(self):
		with h5py.File(f'database/BATH{self.Nb}.h5', 'r') as file:
			self.HGroup = np.array(list(file.keys()))
			if( 0 < self.Hop < 5 ):
				hop = np.vectorize(lambda x: x.startswith(f'N{self.Hop}'))(self.HGroup)
				self.HGroup = self.HGroup[hop]
			elif( self.Hop > 4 ):
				print("args.hop erorr: # of nearest neighbor is wrong")
				sys.exit(1);

			omg, omgn = file['omega']
			if 'omega' in self.HGroup:
				self.HGroup = np.delete(self.HGroup, len(self.HGroup)-1 )
			self.train_list, self.test_list = train_test_split( range(len(self.HGroup)), test_size=0.2, random_state=self.seed )

			Diwn = []; Dw = []; bath = []; chi = []; Gcoeff = []; Sum_rule = []; bath_true = []; #Giwn = []; Dcoeff = [];
			for n,H in enumerate(self.HGroup):
				bath_load = file[f'{H}/bath{self.Nb}/P0'][1]
				bath2 = file[f'{H}/bath{self.Nb}/P0'][1]
				bath_true.append(bath2)
				bath_load[:self.Nb] /= np.sqrt( sum( bath_load[:self.Nb]**2 ) )
				bath.append(bath_load)	
				bath[n][:self.Nb] = np.abs(bath[n][:self.Nb])

				chi.append(file[f'{H}/bath{self.Nb}/P0'].attrs['chi'][0])
				Data = file[H]
					
				Diwn.append( Data['Diwn'][0] + 1j*Data['Diwn'][1] );	#Giwn.append( Data['Giwn'][0] + 1j*Data['Giwn'][1] )
					
				Dw.append( Data['Dw'][0] + 1j*Data['Dw'][1] )	
					
				Gcoeff.append( Data['Gcoeff'][0][:2*self.Lmax] )
		#		Dcoeff.append( Data['Dcoeff'][0][:2*self.Lmax] )

			if( self.Lmax != 0):
				Gcoeff = normalize(Gcoeff, norm="l1")	#normailzation
		
		#	for i in range(self.Lmax):
		#		Gcoeff[:,i] = (Gcoeff[:,i] - Gcoeff[:,i].mean()) / Gcoeff[:,i].std()	# standardization	

		#	혼성함수 자체를 feature로 하는 경우도 해봐야  
		return len(Dw[0]), np.stack(Diwn), np.stack(Dw), np.array(bath), np.array(bath_true), np.array(chi), 1j*omgn, omg, np.array(Gcoeff)#, np.array(Dcoeff), np.stack(Giwn)

	def Pred_opt(self):
		file_name = f'{self.model_path}/BATH{self.Nb}_opt.h5'
		opt_bath = []; opt_chi = [];	
		with h5py.File(file_name,'r') as file:
			for i, group in enumerate(self.HGroup[self.test_list]):
				true_index = np.where(self.HGroup == group)[0][0]
				print(self.HGroup[true_index]);

				Group = file[f'{group}']
				pred_bath = Group['pred_bath']
				pred_chi = pred_bath.attrs['pred_chi']
				opt_bath.append( Group['opt_bath'][0] )
				opt_chi.append( Group['opt_bath'].attrs['opt_chi'][0])
				return opt_bath, opt_chi
				sys.exit(1);
			#	print(pred_bath[0])
			#	print(opt_bath[0]);
			#	print(self.bath_true[true_index])
				print(f'opt: {opt_chi[0]}, true: {self.chi[true_index]} , diff: {abs(opt_chi[0] - self.chi[true_index])}');
		'''			
				fig = plt.figure()
				plt.scatter(pred_bath[0][self.Nb:], pred_bath[0][:self.Nb], s=40, label='pred');
				plt.scatter(opt_bath[0][self.Nb:], opt_bath[0][:self.Nb], s=40, label='opt');
				plt.scatter(self.bath_true[true_index][self.Nb:], np.abs(self.bath_true[true_index][:self.Nb]), s=40, label='true');
			#	plt.title(f'$\\chi_{{opt}}$ - $\\chi_{{true}}$ = {abs(opt_chi[0] - self.chi[self.test_list][i]):.2e}', fontsize = 15)
				plt.title(f'{self.HGroup[true_index]}', fontsize=15)
				plt.xlabel("$\\epsilon_{l}$", fontsize=18)
				plt.ylabel("$V_{l}$", fontsize=18)
				plt.legend(framealpha=0.3)
				plt.show()
			#	fig.savefig(f'{self.model_path}/Nb{self.Nb}_{self.HGroup[true_index]}.png', transparent=True)
				sys.exit(1)
		'''		

	def Chi(self, bath):
		Dw_r = self.Dw[self.test_list]
			 
		#Renolmailzation of V
		for i in range(len(self.test_list)):
#			y = lambda k: np.imag(Dw_r)[i][np.searchsorted(self.w,k).tolist()]
#			Sum, err = integrate.quad( y, self.w[0], self.w[-1] )
#			bath_p[i][:self.Nb] *= np.sqrt(-Sum/np.pi)
			bath[:self.Nb] *= np.sqrt( -integrate.simpson( np.imag(Dw_r[i]), x=self.w ) /np.pi )
		
		Diwn = []
		for i in range(len(self.test_list)):
			Diwn.append([sum(bath[i][l]**2 / (self.iwn[j] - bath[i][self.Nb + l]) for l in range(self.Nb)) for j in range(self.lim)])
		Diwn = np.array(Diwn)
	
		chi = np.sum( np.abs(Diwn - self.Diwn[self.test_list][:, :self.lim])**2, axis=1) / self.lim
#		print(chi[0],'\n',self.chi[self.test_list][0],'\n')
	
		return chi, Diwn

	def Plot_hyb(self):
		fig, ax = plt.subplots(1,1, figsize=(12.8,9.6))

		ax.plot(np.imag(self.iwn)[:self.lim], np.imag(self.Diwn[self.test_list][0])[:self.lim], label='true')
		ax.scatter(np.imag(self.iwn)[:self.lim], np.imag(self.Diwn_p[0])[:self.lim], label=f'${self.Mtype}^p$$_{{L={self.Lmax}}}$')
		ax.scatter(np.imag(self.iwn)[:self.lim], np.imag(self.opt_Diwn[0])[:self.lim], label=f'${self.Mtype}^o$$_{{L={self.Lmax}}}$')
		plt.legend()
		plt.show()

	def Linear_Regression(self):
		X_train = self.Gcoeff[self.train_list]; y_train = self.bath[self.train_list];
		
		model_path = f'{self.model_path}/{self.model_name}'
		if os.path.exists(model_path):
			model = load(model_path)
		else:
			model = LinearRegression()
			model.fit(X_train, y_train)
			dump(model, model_path)  # model save to .joblib
		return model

	def Random_Forest(self):
		X_train = self.Gcoeff[self.train_list]; y_train = self.bath[self.train_list];

#		prams = { 'n_estimators': [200, 2000],
#			'max_depth':[None]
#		}
#		grid_search = GridSearchCV(estimator=RandomForestRegressor(), param_grid=prams, scoring='r2', cv=2)	
#		grid_search.fit(X[self.train_list], y[self.train_list])
		# 최적의 하이퍼파라미터와 예측 정확도 출력
#		print('최적 하이퍼파라미터:', grid_search.best_params_)
#		print('최적 예측 정확도: {0:.4f}'.format(-grid_search.best_score_))

		model_path = f'{self.model_path}/{self.model_name}'
		if os.path.exists(model_path):
			model = load(model_path)
		else:
			model = RandomForestRegressor(n_estimators=1000, random_state=self.seed)
			model.fit(X_train, y_train)
			#dump(model, model_path)  # model save to .joblib
		return model
	
	def Save_type(self, file_name, open_type):
		with h5py.File(file_name, open_type) as file:
			with h5py.File(f'database/BATH{self.Nb}.h5', 'r') as temp:
				dataset = file.create_dataset('omega', data = temp['omega'])
				for i, group in enumerate(self.HGroup[self.test_list]):
					Group = file.create_group(f'{group}')
					reshaped_test_p = np.reshape(self.test_p[i], (1,2*self.Nb))
					dataset = Group.create_dataset('pred_bath', data = reshaped_test_p)
					dataset.attrs['pred_chi'] = self.chi_p[i]
					Diwn = Group.create_dataset('Diwn', data = temp[f'{group}']['Diwn'])
	
	def Save_pred(self):
		file_name = f'{self.model_path}/BATH{self.Nb}_opt.h5'
		if os.path.exists(file_name):
			#self.Save_type(file_name, 'a')	
			print("already exist opt file")
			return; 
		else:
			self.Save_type(file_name, 'w')	

	def Predict(self):
		train_p, test_p = self.model.predict(self.Gcoeff[self.train_list]), self.model.predict(self.Gcoeff[self.test_list])
		return train_p, test_p

	def Plot_sum(self):
		fig, ax = plt.subplots(2,1, figsize=(12.8,9.6))

		ax[0] = plt.subplot(2,1,1); ax[1] = plt.subplot(2,1,2);
		ax[0].set_xlim(-1.5,1.5)
		
		x1 = np.where(self.HGroup == self.HGroup[self.test_list][0] )[0][0]
		x2 = np.where(self.HGroup[self.test_list] == self.HGroup[x1] )[0][0]
	
		pred_scal = self.test_p.copy()
		for i in range( len(pred_scal) ):	# all hybridization V scaling to satisfy sum rule
			pred_scal[i, :self.Nb] *= np.sqrt( 1. / sum(pred_scal[i, :self.Nb]**2) )
	
		bath_true = self.bath[x1]
		bath_true[:self.Nb] = np.abs(bath_true[:self.Nb])
		if (self.plot_already_done == False):
			Dw_true = [sum(bath_true[l]**2 / (self.w[j] + 1j * self.eta - bath_true[self.Nb + l]) for l in range(self.Nb)) for j in range(self.Nmax)]
			ax[0].plot(self.w, np.imag(Dw_true), color='black', label='true')
			ax[1].scatter(bath_true[self.Nb:], bath_true[:self.Nb], color='black', label='true')
	
		bath_p = pred_scal[x2]
		Dw_p = [sum(bath_p[l]**2 / (self.w[j] + 1j * self.eta - bath_p[self.Nb + l]) for l in range(self.Nb)) for j in range(self.Nmax)]
		ax[0].plot(self.w, np.imag(Dw_p), label=f'{self.Mtype}$^s$$_{{L={self.Lmax}}}$')
		ax[1].scatter(bath_p[self.Nb:], bath_p[:self.Nb], label=f'${self.Mtype}^s$$_{{L={self.Lmax}}}$')
		print(f'{self.HGroup[x1]}: {r2_score(bath_true, bath_p)}')
		plt.legend(loc='upper right')
		plt.show()

	def R2_score(self, test):
		print(f'Lmax={self.Lmax:2} \t each R2 score')
		V_r2 = []; e_r2 = []; tot_r2 = [];
		for i in range(len(self.test_list)):
			V_r2.append(r2_score(self.bath[self.test_list[i]][:self.Nb], test[i, :self.Nb]))
			e_r2.append(r2_score(self.bath[self.test_list][i, self.Nb:], test[i, self.Nb:]))
			tot_r2.append(r2_score(self.bath[self.test_list][i], test[i]))
			print(f'H: {self.HGroup[self.test_list][i]:12} V:{V_r2[i]:6.4f} ' \
					f'e:{e_r2[i]: 6.4f} ' \
					f'tot:{tot_r2[i]: 6.4f}')
		np.array(V_r2); np.array(e_r2); np.array(tot_r2);
		print(f'mean R2 score> V:{np.mean(V_r2):6.4f} ' \
				f'e:{np.mean(e_r2): 6.4f} ' \
				f'tot:{np.mean(tot_r2): 6.4f}')

parser = argparse.ArgumentParser();
parser.add_argument('--model', '-m', type=str, dest='model')
parser.add_argument('--Nb', '-Nb', type=int, dest='Nb')
parser.add_argument('--seed', '-s', type=int, dest='seed')
parser.add_argument('--hop', '-n', type=int, dest='hop') #학습시킬 nearest 최대치를 결정 (5면 1~4까지 모두 학습)
#parser.add_argument('--Lmax', '-L', type=int, dest='Lmax')
args = parser.parse_args();

try:
#	subprocess.run(f'make pred_opt', shell=True, check=True)
	for j in range(1):
		model = Machine_Learning(Mtype=args.model, Nb=args.Nb, Lmax=5+2*j, Hop=args.hop, random_state=args.seed)
		print(len(model.test_p), model.test_p[0])
		print(len(model.bath_true[model.test_list]), model.bath_true[model.test_list][0])
#		subprocess.run(f'./pred_opt {args.Nb}', shell=True, check=True)
#		model.Plot_sum()
		model.Plot_hyb()
#		model.Test();
#		model.R2_score(test=model.test_p)
except KeyboardInterrupt:
	print("")
print("END")
