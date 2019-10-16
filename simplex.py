import numpy as np
from fractions import Fraction

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
			if(matrix[i][j] == 1): #may be a basis column
				is_basis = True
				for k in range(len(matrix)):
					if(matrix[k][j] != 0 and k != i):
						is_basis = False
				if(is_basis == True):
					index.append(j)
					break
	return index

#returns the solution, according to the basis
def x_solution(A, b, c, n, m):
    #append b to A
    A = np.concatenate((A,b.reshape((-1,1))),axis = 1)
   
    base_list = []
    for i in range(m):
        cont = 0
        p_j = 0
        if(c[i] == 0):
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

##Main functions

#used in case theres no trivial basis in A
def auxiliar(A, b, original_c, certf):
	n = A.shape[0] #number of restrictions
	m = A.shape[1] #number of variables
	v = Fraction(0) #objective value
	
	c_certf = (np.zeros(A.shape[0], dtype='int')) + Fraction()
	eye = (np.eye(n, dtype='int')) + Fraction()
	ones = (np.ones(n, dtype='int')) + Fraction()

	A = np.append(A, eye, axis=1)

	c = (np.zeros(m, dtype='int')) + Fraction()
	c = np.append(c, ones)

	#transforming the auxiliary program to the canonical form
	for i in range(A.shape[0]):
		v += -b[i]
		c += -A[i][:]
		c_certf += -certf[i][:]

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
		A = A[:,:m]
		basis_check = {}
		basis_row = {}
		one_found = {}
		
		#finding the basis
		for i in range(n):
		    for j in range(m):
		        if j not in basis_check:
		            basis_check[j] = 0
		            one_found[j] =  False
		        if abs(A[i][j]) == Fraction(0):
		            pass
		        elif A[i][j] == Fraction(1) and not one_found.get(j):
		            basis_row[j] = i
		            basis_check[j] -= 1
		            one_found[j] = True    
		        else:
		            basis_check[j] += abs(A[i][j])

		#canonical basis
		for i in range(m):
			if basis_check[i] == -1:
				r = basis_row.get(i)
				d = original_c[i]
		        
				v -= d * b[r]
				c_certf -= d * certf[r,:]
				original_c -= d * A[r,:]

		return "otima", A, b, original_c, v, c_certf, certf
	
	else:
		return "inviavel", A[:,:m], b, c, v, c_certf, certf

def simplex(A, b, c):
	c_certf = (np.zeros(A.shape[0], dtype='int')) + Fraction()
	certf = (np.eye(A.shape[0], dtype='int')) + Fraction()
	status = "otima"
	v = Fraction(0)
	
	A, b, negative, certf = negative_b(A, b, certf)
	if negative:
		status, A, b, c, v, c_certf, certf = auxiliar(A, b, c, certf)
	
	if status == "inviavel":
		return status, "no solution", c_certf, v

	#if all c[j] >= 0, stop.
	while verify_cost_funct(c):
		#getting the index of the lowest value
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

		#performing Gauss operations on the lines relative to the pivot
		for j in range(A.shape[0]):
			if j != i:
				certf[j][:] += (-A[j][k] * certf[i][:])
				b[j] += (-A[j][k] * b[i])	
				A[j][:] += (-A[j][k] * A[i][:])

		v += (-c[k] * b[i])
		c_certf += (-c[k] * certf[i][:])
		c += (-c[k] * A[i][:])
	
	solution = x_solution(A, b, c, n, m)

	return status, solution, c_certf, v

##INPUT:
#c = coeffs. of the objective function to be maximized
#A = restriction matrix
#b = right side of the restrictions

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

n, m = A.shape
original_c = c.copy() * -1
A, b, c = standard_form(A, b, c)

status, solution, certificado, objective = simplex(A, b, c)

##OUTPUT:
#STATUS (Optimal, Unbounded, Infeasible)
#X SOLUTION
#CERTIFICATE

print(status)

if status == 'otima':
	print(float(objective))
	solution = [float(x) for x in solution]
	certificado_list = certificado.tolist()
	certificado_list = [round(float(x), 7) for x in certificado_list]
	
	print(solution)
	print(certificado_list[:m])

elif status == 'ilimitada':
	solution_list = solution.tolist()
	solution_list = [round(float(x), 7) for x in solution_list]
	
	certificado_list = certificado.tolist()
	certificado_list = [round(float(x), 7) for x in certificado_list]
	
	print(solution_list[:m])
	print(certificado_list[:m])

else:
	certificado_list = certificado.tolist()
	certificado_list = [round(float(x), 7) for x in certificado_list]
	print(certificado_list)