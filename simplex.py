##reading input
#obj = função objetivo a ser maximizada
#restr = restricoes a que a funçao objetivo esta sujeita
#b = lado direito das equações de restrição

restr = []
b = []

n, m = map(int, input().split())
obj = list(map(int, input().split()))

for i in range(n):
	eq = list(map(int, input().split()))
	restr.append(eq[:m])
	b.append(eq[m])

## Funções auxiliares
def create_identity(n):
	eye = []
	for i in range(n):
		row = []
		for j in range(n):
			if(i == j): 
				row.append(1)
			else:
				row.append(0) 
		eye.append(row)
	return eye

def verifica_nao_negativo(L):
	for i in range(L):
		if L[i] < 0:
			return false
	return true

## PL auxiliar -> retorna certificado de inviabilidade se v.o < 0
##				  retorna uma base viável se v.o = 0

def pl_auxiliar(n, m, restr, b):
	func_obj = []
	vo = 0
	eye = create_identity(n)
	for i in range(m):
		func_obj.append(0)
	for i in range(n):
		func_obj.append(1)
		for j in range(n):
			restr[i].append(eye[i][j])

pl_auxiliar(n, m, restr, b)