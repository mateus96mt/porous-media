nx = 80
ny = 80

lx = 27
ly = 27

dx = 40
dy = 20

e1 = 1 + dx + dy*nx
e2 = e1 + lx

p1 = 1 + dx + dy*(nx+1)
p2 = p1 + 1
p3 = p1 + nx + 2
p4 = p1 + nx + 1

print("e1: ", e1)
print("e2: ", e2)

print("\np1: ", p1)
print("p2: ", p2)
print("p3: ", p3)
print("p4: ", p4)

print("\nlx: ", lx)
print("ly: ", ly)
