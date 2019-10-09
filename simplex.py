import numpy as np

##reading input
#obj = função objetivo a ser maximizada
#restr = restricoes a que a funçao objetivo esta sujeita
#b = lado direito das equações de restrição

A = []
n, m = map(int, input().split())
c = list(map(int, input().split()))

for i in range(n):
	eq = list(map(int, input().split()))
	A.append(eq)
# b = A[:][m]
# A = A[:][:m]
A = np.array(A)
b = A[:,-1]
A = A[:,:-1]
c = np.array(c) * -1

# print(c)
# print(A)
# print(b)

##Auxiliary functions
def verify_cost_funct(c):
	#simplex (tableau) stopping condition
	if np.min(c) >= 0:
		return False
	else:
		return True

#returns the indexes of the basic variables 
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

#Pre-processing functions
def standard_form(A, b, c):
	n = A.shape[0]
	m = A.shape[1]
	eye = np.eye(n)
	zeros = np.zeros(n)
	
	c = np.append(c, zeros)
	A = np.append(A, eye, axis=1)

	return (A, b, c)

def negative_b(A, b):
	for i in range(len(b)):
		if b[i] < 0:
			A[i][:] = -1 * A[i][:]
	return A, b

#def simplex(N, B, A, b, c):
def simplex(A, b, c):
	#if all c[j] >= 0, stop.
	while verify_cost_funct(c):
		#passou aqui quer dizer que existe pelo menos um negativo, pegar o menor
		k = np.argmin(c)

		#if all variables on the column <= 0, LP is unbounded
		if(np.max(A[:,k]) <= 0):
			print('ilimitada')
			return 'ilimitada'

		#pivot minimizes b[i]/A[i][k], where k is the index of the entering variable
		#pivot index = A[i][k]
		min_val = np.inf
		min_indx = np.NINF
		for i in range(A.shape[0]):
			if A[i][k] > 0:
				if abs(b[i]/A[i][k]) < min_val:
					min_val = abs(b[i]/A[i][k])
					min_indx = i
		i = min_indx

		#set pivot to 1
		A[i][:] = A[i][:] / A[i][k]
		b[i] = b[i] / A[i][k]

		#LINHA DO PIVOT = A[i][:]
		#SUBTRAIR DAS OUTRAS LINHAS, A LINHA DO PIVOT
		for j in range(A.shape[0]):
			if j != i:
				A[j][:] += (-A[j][k] * A[i][:])
				b[j] += (-A[j][k] * b[i])
		c += (-c[k] * A[i][:])


A, b, c = standard_form(A, b, c)
A, b = negative_b(A, b)
simplex(A, b, c)