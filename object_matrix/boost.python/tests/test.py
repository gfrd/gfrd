from numpy import ndarray
from object_matrix import *
import math
from random import shuffle
from decimal import Decimal, Context, ROUND_UP
from sys import getrefcount

def float_cmp(lhs, rhs, tolerance):
    if lhs - rhs < tolerance or rhs - lhs < tolerance or lhs == rhs:
        return 0
    if lhs < rhs:
        return -1
    else:
        return 1

c = ObjectContainer(1.0, 10)
assert c.matrix_size == 10
assert c.cell_size == 0.1
assert len(c) == 0
c[0] = Sphere((0.5, 0.3, 0.2), 0.08)
assert len(c) == 1
assert float_cmp(c[0].x, 0.5, 1e-8) == 0 and \
       float_cmp(c[0].y, 0.3, 1e-8) == 0 and \
       float_cmp(c[0].z, 0.2, 1e-8) == 0 and \
       float_cmp(c[0].radius, 0.08, 1e-8) == 0 and \
       c[0].id == 0
assert c[1] == None
assert isinstance(c[0], SphereRef)

tmp = ndarray(shape = (3, ))
tmp[0] = 0.05
tmp[1] = 0.3
tmp[2] = 0.9
c[1] = Sphere(tmp, 0.05)
assert len(c) == 2
assert float_cmp(c[1].x, 0.05, 1e-8) == 0 and \
       float_cmp(c[1].y, 0.3, 1e-8) == 0 and \
       float_cmp(c[1].z, 0.9, 1e-8) == 0 and \
       float_cmp(c[1].radius, 0.03, 1e-8) == 0 and \
       c[1].id == 1
assert isinstance(c[1], SphereRef)
assert c[2] == None

c[[2]] = Sphere((0.4, 0.3, 0.3), 0.08)
assert len(c) == 3
assert float_cmp(c[[2]].x, 0.4, 1e-8) == 0 and \
       float_cmp(c[[2]].y, 0.3, 1e-8) == 0 and \
       float_cmp(c[[2]].z, 0.3, 1e-8) == 0 and \
       float_cmp(c[[2]].radius, 0.08, 1e-8) == 0 and \
       c[[2]].id == [2] 
assert isinstance(c[[2]], SphereRef)
assert c[3] == None

#for i in c.iterneighbors(Sphere((0.45, 0.23, 0.13), 0.09)):
#    print i
#for i in c.iterneighbors_cyclic(Sphere((0.9, 0.23, 0.0), 0.1)):
#    print i
a = c.neighbors_array(Sphere((0.45, 0.23, 0.13), 0.09))
print a[1][0]
assert getrefcount(a[0]) ==  2
assert getrefcount(a[0][0]) == 2
assert len(a[0]) == 1
assert len(a[1]) == 1
assert a[0][0].id == 0
assert float_cmp(a[0][0].x, c[0].x, 1e-8) == 0 and \
       float_cmp(a[0][0].y, c[0].y, 1e-8) == 0 and \
       float_cmp(a[0][0].z, c[0].z, 1e-8) == 0 and \
       float_cmp(a[0][0].radius, c[0].radius, 1e-8) == 0 and \
       a[0][0].id == c[0].id
a = c.all_neighbors_array((0.45, 0.23, 0.2))
assert float_cmp(
    a[1][0],
    math.sqrt((0.5 - 0.45) ** 2 + (0.3 - 0.23) ** 2 + (0.2 - 0.2) ** 2),
    1e-8) == 0

assert len(c) == 3
del c[0]
assert c[0] == None
assert len(c) == 2

c = ObjectContainer(10, 10)
for i in xrange(0, 1000):
    c[str(i)] = Sphere((0, 0, 0), 0)

l = range(0, 1000)
shuffle(l)
for i, k in enumerate(l):
    assert len(c) == 1000 - i
    del c[str(k)]

c = ObjectContainer(10, 10)
os = 0.25
ss = float(
    Decimal(
        str((math.sqrt(c.cell_size ** 2 + c.cell_size * 2) - os * 2) / 2)
        ).normalize(Context(prec = 2, rounding = ROUND_UP))
    )
for i in xrange(0, 10):
    for j in xrange(0, 10):
        for k in xrange(0, 10):
            key = k + j * 10 + i * 100
            c[key] = Sphere((i + 0.5, j + 0.5, k + 0.5), os)

pat = (8, 4, 1, 4)
p = 1
for i in xrange(0, 35):
    assert pat[i % len(pat)] == len(c.neighbors_array(Sphere((p, p, p), ss))[0])
    p += 0.25
assert pat[37 % len(pat)] != len(c.neighbors_array(Sphere((p, p, p), ss))[0])

p = 1
for i in xrange(0, 40):
    assert pat[i % len(pat)] == len(c.neighbors_array_cyclic(Sphere((p, p, p), ss))[0])
    p += 0.25



