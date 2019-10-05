import numpy as np

##reading input
#obj = função objetivo a ser maximizada
#restr = restricoes a que a funçao objetivo esta sujeita
#b = lado direito das equações de restrição

restr = []
n, m = map(int, input().split())
obj = list(map(int, input().split()))

for i in range(n):
	eq = list(map(int, input().split()))
	restr.append(eq)

## Funções auxiliares
def verify_cost_funct(N, c):
	for j in N:
		if(c[j] > 0):
			return True
	return False

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

#INPUT: l = index of the leaving variable x_l
#		e = index of the entering variable x_e
#OUTPUT: the tuple representing the new slack form
def pivot(N, B, A, b, c, v, l, e):
	#Coeff. of the eq. for the new basic variable x_e
	b_e = b[l]/A[l][e]
	for j in range(N):
		if(N[j] != e):
			A[e][j] = A[l][j]/A[l][e]
	A[e][l] = 1/A[l][e]

	#Coeff. of the remaining constraints
	for i in range(B):
		if(B[i] != l):
			b[i] -= A[i][e]*A[e][j]
			for j in range(N):
				if(N[j] != e):
					A[i][j] -= A[i][e]*A[e][j]
			A[i][l] = -A[i][e]*A[e][j]

	#Objective function
	v += v + c[e]*b[e]
	for j in range(N):
		if(N[j] != e):
			c[j] -= c[e]*A[e][j]
	c[l] = -c[e]*A[e][l]
	
	index_e_in_n, = np.where(N == e)
	index_l_in_b, = np.where(B == l)
	N[index_e_in_n] = l
	B[index_l_in_b] = e
	
	return (N, B, A, b, c, v)

#INPUT: A = restr. matrix, b = value of each equation, c = objective function coeff.
#OUTPUT (N, B, A, b, c, v): N = non-basic indexes, B = basic indexes,				
#							A = restriction matrix, b = value of each equation		
#							c = objective function coeff., v = objective value 		

def initialize_simplex(A, b, c):

	return (N, B, A, b, c, v)


def simplex(A, b, c):
	(N, B, A, b, c, v) = initialize_simplex(A, b, c)
	delta = np.zeros(A.shape[1])
	#stopping condition: if all the variables not in the base <= 0
	while verify_cost_funct(N, c):
		max_value = np.NINF
		e = -1
		#select the entering variable
		for i in N:
			if(c[i] > max_value):
				max_value = c[i]
				e = i
		
		for i in B:
			if A[i][e] > 0:
				delta[i] = b[i]/A[i][e]
			else
				delta[i] = np.inf
		
		min_value = np.inf
		l = -1
		#select leaving variable
		for i in B:
			if(delta[i] < min_value):
				min_value = delta[i]
				l = i
		if(delta[l] == np.inf):
			return "ilimitada" #calcular certificado de ilimitada
		else:
			(N, B, A, b, c, v) = pivot(N, B, A, b, c, v, l, e)