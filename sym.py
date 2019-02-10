from sympy import *
z,t,U,e = symbols("z t U e")
H = Matrix(
    [   
        [  U, -t,  t,  0],
        [ -t,  0,  0, -t],
        [  t,  0,  0,  t],
        [  0, -t,  t,  U]
    ])
    
z = diag(z, z, z, z)
G = (z-H).inv()

phi_gs = Matrix([1,e,-e,1])

print(G)
#print(expand(G.trace()))
#print(H.eigenvals())
#print(H.eigenvects())