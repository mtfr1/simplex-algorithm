import numpy as np

##reading input
#obj = função objetivo a ser maximizada
#restr = restricoes a que a funçao objetivo esta sujeita
#b = lado direito das equações de restrição

#FALTA MULTIPLICAR C POR -1
restr = []
n, m = map(int, input().split())
obj = list(map(int, input().split()))

for i in range(n):
	eq = list(map(int, input().split()))
	restr.append(eq)

##Auxiliary functions
def verify_cost_funct(c):
	#simplex (tableau) stopping condition
	if np.min(c) >= 0:
		return False
	else:
		return True

def standard_form(A, b, c):
	n = A.shape[0]
	m = A.shape[1]
	eye = np.eye(n)
	zeros = np.zeros(n)
	
	c = np.append(c, zeros)
	A = np.append(A, eye, axis=1)

	return (A, b, c)

#encontra os indices das colunas que pertencem a base 
def basis_index(matrix):
	index = []
	for i in range(len(matrix)):
		for j in range(len(matrix[0])):
			if(matrix[i][j] == 1): #possivel coluna pertencente a base
				is_basis = True
				for k in range(len(matrix)):
					if(matrix[k][j] != 0 and k != i):
						is_basis = False
				if(is_basis == True):
					index.append(j)
					break
	return index


def simplex(N, B, A, b, c):
	#tableau stopping condition: if all c[j] >= 0, stop.
	while verify_cost_funct(c):
		max_c = np.NINF
		max_indx = -1
		#select entering variable, index = k
		for k in N:
			if c[k] > max_c:
				max_indx = k
				max_c = c[k]
		k = max_indx

		#if all variables on the column <= 0, LP is unbounded
		if(np.max(A[:][k]) <= 0):
			return 'ilimitada'
		
		#pivot minimizes b[i]/A[i][k], where k is the index of the entering variable
		#pivot index = A[i][k]
		for i in range(len(b)):
			min_val = np.inf
			min_indx = -1
			if A[i][k] > 0:
				if b[i]/A[i][k] < min_val:
					min_val = b[i]/A[i][k]
					min_indx = i
		i = min_indx

		#leaving variable index = leaving_index
		# B - leaving_index + k
		# N - k + leaving_index
		leaving_index = -1
		for e in range(A.shape[1]):
			if(A[i][e] == 1 and e != k):
				is_basis = True
				for l in range(A.shape[0]):
					if(A[l][e] != 0 and l != i):
						is_basis = False
				if(is_basis == True):
					leaving_index = l
					break

		#set pivot to 1
		A[i][:] = A[i][:] / A[i][k]