import numpy




if __name__ == "__main__":
    p = 10
    q = 3
    c = 6.4

    a_ell = 0.5*(p+q)
    b_ell = numpy.sqrt(a_ell**2 - c**2)
    print("Ellipse   a,b,c: ", a_ell, b_ell, c )
    theta_ell = numpy.arcsin(b_ell / numpy.sqrt(p*q))
    print("Ellipse theta: ", theta_ell * 180 /numpy.pi  )
    theta_ell2 = 0.5 * numpy.arccos(-(p ** 2 + q ** 2 - 4 * c ** 2) / (2 * p * q))
    print("Ellipse theta2: ", theta_ell2 * 180 /numpy.pi  )

    a_hyp = 0.5*(p-q)
    b_hyp = numpy.sqrt(c**2 - a_hyp**2)

    print("Hyperbola a,b,c: ", a_hyp, b_hyp, c)


    theta_hyp = 0.5 * numpy.arccos( (p**2 + q**2 - 4 * c**2)  / (2 * p * q))


    print("Hyperbola theta: ", theta_hyp * 180 / numpy.pi)