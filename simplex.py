import numpy as np
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

A = np.array(A, dtype=np.float64)
b = A[:,-1]
A = A[:,:-1]
c = np.array(c, dtype=np.float64) * -1

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

##Pre-processing functions
def standard_form(A, b, c):
	n = A.shape[0]
	m = A.shape[1]
	eye = np.eye(n)
	zeros = np.zeros(n)
	
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
	v = 0 #objective value
	
	c_certf = np.zeros(n)
	eye = np.eye(n)
	ones = np.ones(n)
	
	A = np.append(A, eye, axis=1)
	c = np.zeros(m)
	c = np.append(c, ones)
	
	#transforming the auxiliary program to the canonical form
	for i in range(A.shape[1]-m):
		# c += -c[i+m]*A[i][:]
		# v += -c[i+m]*b[i]
		# c_certf += -c[i+m]*certf[i][:]
		c += -A[i][:]
		v += -b[i]
		c_certf += -certf[i][:]
	
	#now, execute the same steps as the simplex
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
		A[i][:] = A[i][:] / A[i][k]
		b[i] = b[i] / A[i][k]

		for j in range(A.shape[0]):
			if j != i:
				certf[j][:] += (-A[j][k] * certf[i][:])
				A[j][:] += (-A[j][k] * A[i][:])
				b[j] += (-A[j][k] * b[i])

		c_certf += (-c[k] * certf[i][:])
		c += (-c[k] * A[i][:])
		v += (-c[k] * b[i])
	
	if(v == 0):
		return "otima", A[:,:m], b, c_certf, certf
	else:
		return "inviavel", A[:,:m], b, c_certf, certf

def simplex(A, b, c):
	c_certf = np.zeros(A.shape[0])
	certf = np.eye(A.shape[0])
	status = "otima"

	A, b, negative, certf = negative_b(A, b, certf)
	if negative:
		status, A, b, c_certf, certf = auxiliar(A, b, certf)
	
	if status == "inviavel":
		return status, "no solution", c_certf

	#if all c[j] >= 0, stop.
	while verify_cost_funct(c):
		#passou aqui quer dizer que existe pelo menos um negativo, pegar o menor
		k = np.argmin(c)

		#if all variables on the column <= 0, LP is unbounded
		if(np.max(A[:,k]) <= 0):
			status = 'ilimitada'
			solution = np.zeros(A.shape[1])
			c_certf = np.zeros(A.shape[1])
			c_certf[k] = 1
			
			basis = basis_index(A)
			for j in basis:
				for i in range(A.shape[0]):
					if A[i][j] == 1:
						c_certf[j] = -A[i][k]
						solution[j] = b[i]

			return status, solution, c_certf

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
		A[i][:] = A[i][:] / A[i][k]
		b[i] = b[i] / A[i][k]

		#LINHA DO PIVOT = A[i][:]
		#SUBTRAIR DAS OUTRAS LINHAS, A LINHA DO PIVOT
		for j in range(A.shape[0]):
			if j != i:
				certf[j][:] += (-A[j][k] * certf[i][:])
				A[j][:] += (-A[j][k] * A[i][:])
				b[j] += (-A[j][k] * b[i])	
		c_certf += (-c[k] * certf[i][:])
		c += (-c[k] * A[i][:])
	
	B = basis_index(A)
	solution = np.zeros(A.shape[0])
	for i in range(A.shape[0]):
		if i in B:
			solution[i] = b[i]
		else:
			solution[i] = 0

	return status, solution, c_certf


original_c = c.copy() * -1
A, b, c = standard_form(A, b, c)

status, solution, certificado = simplex(A, b, c)

print(status)
if status == 'otima':
	v = np.sum(original_c * solution[:m])
	print(v)
	print(solution[:m])
	print(certificado[:m])

elif status == 'ilimitada':
	print(solution[:m])
	print(certificado[:m])