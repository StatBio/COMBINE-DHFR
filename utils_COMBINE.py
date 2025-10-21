####################### MODULES #######################
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
from config import CODE, AA_2_NUM, POS_ALIGN_Thompson, POS_ALIGN_Kalmer, SEQ_ECOLI_Kalmer, SEQ_ECOLI_Thompson

#######################################################

####################### Statistics ####################

def Stats_Pos_bis(align, Pos,ref_pos=None, maxp=5,perc=True,MSA_name='Thompson'):
	if MSA_name == 'Thompson':
		POS_ALIGN = POS_ALIGN_Thompson
	elif MSA_name == 'Kalmer':
		POS_ALIGN = POS_ALIGN_Kalmer
	else: 
		POS_ALIGN = np.arange(1,align.shape[1]+1)
	
	list_aa = align[:, np.where(POS_ALIGN == Pos)[0]]

	if ref_pos is not None:
		list_aa_ref = align[:, np.where(POS_ALIGN == ref_pos[1])[0]]
		ind_ref = (list_aa_ref==AA_2_NUM[ref_pos[0]])
		list_aa_withref = list_aa[ind_ref]
		unique, counts = np.unique(list_aa_withref, return_counts=True)
	else:
		unique, counts = np.unique(list_aa, return_counts=True)

	sorted_data = sorted(zip(counts, unique), reverse=True)
	counts, unique = zip(*sorted_data)
	if len(counts) > maxp:
		top_counts = list(counts[:maxp])
		top_labels = [CODE[aa] for aa in unique[:maxp]]
		others_count = sum(counts[maxp:])
		top_counts.append(others_count)
		top_labels.append('Others')
	else:
		top_counts = list(counts)
		top_labels = [CODE[i] for i in unique]

	# Trier les données par ordre décroissant
	sorted_data = sorted(zip(top_counts, top_labels), reverse=True)
	top_counts, top_labels = zip(*sorted_data)
	top_counts,top_labels = list(top_counts),list(top_labels)

	percentages = [100 * top_counts[i] / sum(top_counts) for i in range(len(top_counts)) if top_labels[i]!='Others']
	top_counts = [top_counts[i] for i in range(len(top_counts)) if top_labels[i]!='Others']
	top_labels = [top_labels[i] for i in range(len(top_labels)) if top_labels[i]!='Others']

	plt.figure(figsize=(6, 4))
	bars = plt.barh(range(len(top_counts)), top_counts, color='black')

	plt.ylabel("Amino Acids", fontsize=18)
	plt.yticks(range(len(top_labels)), top_labels, fontsize=16)

	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.gca().spines['bottom'].set_visible(False)

	plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

	for i, (count, percentage) in enumerate(zip(top_counts, percentages)):
		if perc:
			plt.text(count, i, f'  {percentage:.1f}%', va='center',fontsize=16)
		else: plt.text(count, i, f'  {count}', va='center',fontsize=16)

	plt.gca().invert_yaxis() 
	plt.tight_layout()

	plt.show()

#######################################################


####################### COMBINE #######################


def Propose_Mutation_DDE2(file,path,Nb_mut,Mutations,Np=10,MSA_name='Thompson',SEQ_REF = None):
	output = np.load(path+file+'.npy',allow_pickle=True)[()]
	J,h = Jw(output['W_all'][0],output['options1']['q'])
	J,h = Zero_Sum_Gauge(J,h)

	if MSA_name == 'Thompson':
		POS_ALIGN = POS_ALIGN_Thompson
		if SEQ_REF is None:
			SEQ_REF = SEQ_ECOLI_Thompson
	elif MSA_name == 'Kalmer':
		POS_ALIGN = POS_ALIGN_Kalmer
		if SEQ_REF is None:
			SEQ_REF = SEQ_ECOLI_Kalmer
	else: 
		POS_ALIGN = np.arange(1,h.shape[0]+1)
		if SEQ_REF is None:
			raise ValueError("A reference sequence (SEQ_REF) must be provided for prediction.")
	
	if len(SEQ_REF) != h.shape[0]:
		raise ValueError("The alignment length must match the length of the reference sequence.")

	seq_tryp = np.array([AA_2_NUM[SEQ_REF[i]] for i in range(len(SEQ_REF))])

	Ind_proposed = []
	for Mut in Mutations:
		seq_mut = np.copy(seq_tryp)
		seq_mut[np.where(POS_ALIGN==Mut[1])[0][0]] = AA_2_NUM[Mut[0]]

		Ind_proposed.append((AA_2_NUM[Mut[0]],np.where(POS_ALIGN==Mut[1])[0][0]))

	Mut_pred_WT = Mutational_effect(seq_tryp,h,J)

	seq_ref = np.copy(seq_mut)

	if Nb_mut>0:
		for mut in range(Nb_mut):

			Mut_pred = Mutational_effect(seq_ref,h,J)

			Mut_diff = Mut_pred - Mut_pred_WT

			for ind in Ind_proposed:
				Mut_diff[:,ind[1]] = np.nan

			if mut==0:
				plot_hist(Mut_diff.flatten())

			ind_min = np.unravel_index(np.nanargmin(Mut_diff, axis=None), Mut_diff.shape)
			Ind_proposed.append(ind_min)
			sig_mut = ind_min[0]

			if mut < Np:
				print('Mutation '+str(mut + 2)+': '+str(CODE[seq_ref[ind_min[1]]])+str(POS_ALIGN[ind_min[1]])+str(CODE[sig_mut]),Mut_diff[ind_min],Mut_pred[ind_min],Mut_pred_WT[ind_min])
			seq_ref[ind_min[1]] = sig_mut

	return Ind_proposed

def plot_hist(data,ylimit=(0,40)):
	fig, ax = plt.subplots(figsize=(6, 4))

	ax.hist(data, bins=50, color='grey', edgecolor='black', alpha=0.7)
	ax.set_ylim(ylimit)
	ax.set_xlabel(r'$\Delta \Delta E_{\mathrm{EcDHFR}}^{\mathrm{G121V}}$', fontsize=18)
	ax.set_ylabel('Occurences', fontsize=18)
	ax.grid(axis='y', linestyle='--', alpha=0.7)

	ax.tick_params(axis='both', labelsize=18)
	for spine in ["top", "right"]:
		ax.spines[spine].set_visible(False)

	plt.tight_layout()

	plt.show()

#######################################################

####################### TOOLS #######################

def Zero_Sum_Gauge(J,h):
	"""
	Function to apply a zero-sum gauge transformation to J and h matrices.
	"""
	J_zg = np.copy(J)
	h_zg = np.copy(h)

	h_zg -= np.expand_dims(np.mean(h,axis = 1),axis = 1) 
	h_zg += np.sum(np.mean(J,axis=3)-np.expand_dims(np.mean(J,axis=(2,3)),axis=2),axis=1)

	J_zg -= np.expand_dims(np.mean(J,axis = 2),axis = 2) 
	J_zg -= np.expand_dims(np.mean(J,axis=3),axis =3) 
	J_zg += np.expand_dims(np.mean(J,axis=(2,3)),axis=(2,3))
	return J_zg, h_zg

def Mutational_effect(WT_seq,h,J):
	"""
	Function to compute mutational effects based on the provided wild-type sequence, h, and J values, 
	considering the changes in energy resulting from single amino acid substitutions.
	"""
	L,q = h.shape
	Mut = np.zeros((q,L))
	for k in range(L):
		for b in range(q):
			if b != WT_seq[k]:
				Delta_E = h[k,WT_seq[k]] - h[k,b] - np.sum(J[k,np.arange(L),b,WT_seq]) + np.sum(J[k,np.arange(L),WT_seq[k],WT_seq]) 
				Mut[b,k] = Delta_E
	return Mut

def Jw(W,q):
	L=int(((q*q-2*q)+((2*q-q*q)**2+8*W.shape[0]*q*q)**(1/2))/2/q/q)
	J=np.zeros((L,L,q,q))
	h=np.zeros((L,q))
	x=np.array([[i,j] for i,j in it.combinations(range(L),2)])
	for a in range(q):
		for b in range(q):
			J[x[:,0],x[:,1],a,b]=W[(q**2*((x[:,0])*(2*L-x[:,0]-1)/2+(x[:,1]-x[:,0]-1))+(a)*q+b).astype(int)]
			J[x[:,1],x[:,0],b,a]=J[x[:,0],x[:,1],a,b]
	x=np.array(range(L))
	for a in range(q):
		h[x[:],a]=W[(q**2*L*(L-1)/2+q*x[:]+a).astype(int)]
	return J,h

#######################################################