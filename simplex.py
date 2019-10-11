import numpy as np
from fractions import Fraction
np.set_printoptions(formatter={"float_kind": lambda x: "%g" % x})

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

A = np.array(A)
b = A[:,-1]
A = A[:,:-1]
c = np.array(c) * -1

A = A + Fraction()
b = b + Fraction()
c = c + Fraction()

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
def round_array(c):
	for i in range(len(c)):
		if(abs(0 - c[i]) < 1e-6):
			c[i] = 0
	return c

def printa_bonitim(A):
	for line in A:
		fline = list(map(float, line))
		print(fline)
def find_basis(A, b, c, n, m):
    #Append b to A
    A = np.concatenate((A,b.reshape((-1,1))),axis = 1)
   
    base_list = []
    for i in range(m):
        cont = 0
        p_j = 0;
        if(c[i]==0):
            for j in range(n):
                if (A[j][i] == 1):
                    p_j = j
                    cont += 1
        if(cont != 1):
            base_list.append(0)
        if(cont == 1):
            base_list.append(b[p_j])
    return base_list

##Pre-processing functions
def standard_form(A, b, c):
	n = A.shape[0]
	m = A.shape[1]
	eye = (np.eye(n, dtype='int')) + Fraction()
	zeros = (np.zeros(n, dtype='int')) + Fraction()
	
	c = np.append(c, zeros)
	A = np.append(A, eye, axis=1)

	return (A, b, c)

def negative_b(A, b, certf):
	negative = False
	for i in range(len(b)):
		if b[i] < 0:
			negative = True
			A[i][:] = -1 * A[i][:]
			b[i] = -1 * b[i]
			certf[i][:] = -1 * certf[i][:]
	return A, b, negative, certf


#used in case theres no trivial basis in A
def auxiliar(A, b, certf):
	n = A.shape[0] #number of restrictions
	m = A.shape[1] #number of variables
	v = Fraction(0) #objective value
	
	c_certf = (np.zeros(A.shape[0], dtype='int')) + Fraction()
	eye = (np.eye(n, dtype='int')) + Fraction()
	ones = (np.ones(n, dtype='int')) + Fraction()

	A = np.append(A, eye, axis=1)

	c = (np.zeros(m, dtype='int'))+ Fraction()
	c = np.append(c, ones)

	#transforming the auxiliary program to the canonical form
	for i in range(A.shape[0]):
		v += -b[i]
		c += -A[i][:]
		c_certf += -certf[i][:]
	#CORRETO ATE AQUI
	
	#execute the same steps as the simplex
	while verify_cost_funct(c):
		k = np.argmin(c)

		min_val = np.inf
		min_indx = np.NINF
		for i in range(A.shape[0]):
			if A[i][k] > 0:
				if abs(b[i]/A[i][k]) < min_val:
					min_val = abs(b[i]/A[i][k])
					min_indx = i
		i = min_indx


		certf[i][:] = certf[i][:] / A[i][k]
		b[i] = b[i] / A[i][k]
		A[i][:] = A[i][:] / A[i][k]
		
		
		for l in range(A.shape[0]):
			if l != i:
				certf[l][:] += (-A[l][k] * certf[i][:])
				b[l] += (-A[l][k] * b[i])
				A[l][:] += (-A[l][k] * A[i][:])

		v += (-c[k] * b[i])
		c_certf += (-c[k] * certf[i][:])
		c += (-c[k] * A[i][:])

	if(v == 0):
		return "otima", A[:,:m], b, c_certf, certf
	else:
		return "inviavel", A[:,:m], b, c_certf, certf

def simplex(A, b, c):
	c_certf = (np.zeros(A.shape[0], dtype='int')) + Fraction()
	certf = (np.eye(A.shape[0], dtype='int')) + Fraction()
	status = "otima"
	v = Fraction(0)
	
	A, b, negative, certf = negative_b(A, b, certf)
	if negative:
		status, A, b, c_certf, certf = auxiliar(A, b, certf)
	
	if status == "inviavel":
		return status, "no solution", c_certf, v

	
	#if all c[j] >= 0, stop.
	while verify_cost_funct(c):
		#passou aqui quer dizer que existe pelo menos um negativo, pegar o menor
		k = np.argmin(c)

		#if all variables on the column <= 0, LP is unbounded
		if(np.max(A[:,k]) <= 0):
			status = 'ilimitada'
			solution = (np.zeros(A.shape[1], dtype='int')) + Fraction()
			c_certf = (np.zeros(A.shape[1], dtype='int')) + Fraction()
			c_certf[k] = 1
			
			basis = basis_index(A)
			for j in basis:
				for i in range(A.shape[0]):
					if A[i][j] == 1:
						c_certf[j] = -A[i][k]
						solution[j] = b[i]

			return status, solution, c_certf, v

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
		certf[i][:] = certf[i][:] / A[i][k]
		b[i] = b[i] / A[i][k]
		A[i][:] = A[i][:] / A[i][k]

		#LINHA DO PIVOT = A[i][:]
		#SUBTRAIR DAS OUTRAS LINHAS, A LINHA DO PIVOT
		for j in range(A.shape[0]):
			if j != i:
				certf[j][:] += (-A[j][k] * certf[i][:])
				b[j] += (-A[j][k] * b[i])	
				A[j][:] += (-A[j][k] * A[i][:])

		v += (-c[k] * b[i])
		c_certf += (-c[k] * certf[i][:])
		c += (-c[k] * A[i][:])
	
	# B = basis_index(A)
	# solution = (np.zeros(A.shape[1]-n, dtype='int')) + Fraction()
	
	# # for i in range(len(b)):
	# #  	solution[B[i]] = b[i]
	# #  	# if c[i] == 0 and i < len(b):
	# #  	# 	solution[i] = b[i]
	# #  	# else:
	# #  	# 	solution[i] = 0

	# # for j in range(len(c)):
	# # 	if c[j] == 0:
	# # 		n_zero = 0
	# # 		n_um = 0
	# # 		for i in range(A.shape[0]):
	# # 			if A[i][j] == 1:
	# # 				n_um += 1
	# # 				index = j
	# # 			if A[i][j] == 0:
	# # 				n_zero += 1
	# # 		if n_zero == A.shape[0]-1 and n_um == 1:
	# # 			solution[j] = b[index]
	# # print(solution)
	solution = find_basis(A, b, c, n, m)

	return status, solution, c_certf, v

n, m = A.shape
original_c = c.copy() * -1
A, b, c = standard_form(A, b, c)

status, solution, certificado, objective = simplex(A, b, c)

print(status)

if status == 'otima':
	print(float(objective))
	solution = [float(x) for x in solution]
	
	certificado_list = certificado.tolist()
	certificado_list = [float(x) for x in certificado_list]
	
	print(solution)
	print(certificado_list[:m])

elif status == 'ilimitada':
	solution_list = solution.tolist()
	solution_list = [float(x) for x in solution_list]
	certificado_list = certificado.tolist()
	certificado_list = [float(x) for x in certificado_list]
	
	print(solution_list[:m])
	print(certificado_list[:m])

else:
	certificado_list = certificado.tolist()
	certificado_list = [float(x) for x in certificado_list]
	print(certificado_list)