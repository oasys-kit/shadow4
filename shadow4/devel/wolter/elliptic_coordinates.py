
# see doi: 10.11648/j.mma.20170204.12

import numpy


def ell2cart(q1,q2, c=1.0, verbose=0):
    x = c * numpy.cosh(q1) * numpy.cos(q2)
    y = c * numpy.sinh(q1) * numpy.sin(q2)
    if verbose:
        print("  ")
        print("c: ", c)
        a_hyp, b_hyp = c * numpy.cos(q2), c * numpy.sin(q2)
        c_hyp = numpy.sqrt(a_hyp**2 + b_hyp**2)
        print("hyp a,b,c: ",  a_hyp, b_hyp, c_hyp)
        a_ell, b_ell = c * numpy.cosh(q1), c * numpy.sinh(q1)
        c_ell = numpy.sqrt(a_ell ** 2 + b_ell ** 2)
        print("ell a,b: ", a_ell, b_ell, c_ell)
    return x, y

def cart2ell(x, y, c=1.0, verbose=1):
    B = x**2 + y**2 - c**2
    p = (-B + numpy.sqrt(B**2 + (2 * c * y)**2)) / (2 * c**2)
    q = (-B - numpy.sqrt(B**2 + (2 * c * y)**2)) / (2 * c**2)

    q2 = numpy.arcsin(numpy.sqrt(p))

    if x >= 0 and y >= 0:
        pass
    elif x < 0 and y >= 0:
        q2 = numpy.pi - q2
    elif x <= 0 and y < 0:
        q2 = numpy.pi + q2
    elif x > 0 and y < 0:
        q2 = 2 * numpy.pi - q2

    q1 = 0.5 * numpy.log(1 - 2 * q + 2 * numpy.sqrt(q**2 - q))

    if verbose:
        print("  ")
        print("c: ", c)
        a_hyp, b_hyp = c * numpy.cos(q2), c * numpy.sin(q2)
        c_hyp = numpy.sqrt(a_hyp**2 + b_hyp**2)
        print("hyp a,b,c: ",  a_hyp, b_hyp, c_hyp)
        a_ell, b_ell = c * numpy.cosh(q1), c * numpy.sinh(q1)
        c_ell = numpy.sqrt(a_ell ** 2 + b_ell ** 2)
        print("ell a,b: ", a_ell, b_ell, c_ell)

    return q1, q2



# print(ell2cart(1.0, numpy.pi/2))
# print(ell2cart(1.0, 0.0))
# print(ell2cart(0.0, numpy.pi/2))

# for i in range(10):
#     tmp_ell = (numpy.random.rand() * 5, numpy.random.rand() * 2 * numpy.pi) #TODO: discrepancies for > 5
#     # tmp_ell = (20.0, 300 * numpy.pi / 180)
#     print("  ")
#     print("in (ell): ",tmp_ell)
#     tmp_cart = ell2cart(tmp_ell[0],tmp_ell[1])
#     print("out (cart): ",tmp_cart)
#     tmp_again = cart2ell(tmp_cart[0], tmp_cart[1])
#     print("again (ell): ",tmp_again[0], tmp_again[1])
#     print (numpy.abs(tmp_ell[0]-tmp_again[0]) )
#     print (numpy.abs(tmp_ell[1]-tmp_again[1]) )
#     assert (numpy.abs(tmp_ell[0]-tmp_again[0]) < 1e-10 )
#     assert (numpy.abs(tmp_ell[1]-tmp_again[1]) < 1e-10 )


ssour=10
simag=3
theta_grazing=3e-3

a = 0.5 * (ssour + simag)
b = numpy.sqrt(ssour * simag) * numpy.sin(theta_grazing)
c = numpy.sqrt(a**2 - b**2)

YCEN = (ssour**2 - simag**2) / 4 / c
ZCEN = -b * numpy.sqrt(1 - YCEN**2 / a**2)


print("ssour,simag,c,YCEN,ZCEN", ssour,simag,c,YCEN,ZCEN)
print("a,b,c (ell)", a, b, c)

print("Yc,Zc: ", YCEN, ZCEN)
q1, q2 = cart2ell(YCEN, ZCEN, c=c)
print("q1,q2: ", q1, q2)
Yc, Zc = ell2cart(q1, q2, c=c)
print("Yc,Zc: ", Yc, Zc )
d1 = numpy.sqrt( (Yc + c)**2 + Zc**2) + numpy.sqrt( (Yc - c)**2 + Zc**2)
print(">>>>", d1,  2* a, ssour + simag)

# another point on the ellipse
Yc2, Zc2 = ell2cart(q1, q2 * 1.2, c=c)
print("Yc2,Zc2: ", Yc2, Zc2)
d2 = numpy.sqrt( (Yc2 + c)**2 + Zc2**2) + numpy.sqrt( (Yc2 - c)**2 + Zc2**2)
print(">>>>", d2,  2* a, ssour + simag)

# another point on the hyperbola
Yc3, Zc3 = ell2cart(q1 * 1.2, q2, c=c, verbose=1)
print("Yc3,Zc3: ", Yc3, Zc3)
d3 = numpy.sqrt( (Yc3 + c)**2 + Zc3**2) + numpy.sqrt( (Yc3 - c)**2 + Zc3**2)
print(">>>>", d3,  2* a, ssour + simag)