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
