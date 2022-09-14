import numpy as np

def rotate_and_translate_coefficients(coe_list,R_M,T):
    axx, ayy, azz, axy, ayz, axz, ax, ay, az, a0 = coe_list

    A2 = np.array([[axx,axy/2,axz/2],
                   [axy/2,ayy,ayz/2],
                   [axz/2,ayz/2,azz]])
    A1 = np.array([ax,ay,az])
    A0 = a0

    B2 = np.dot(R_M, np.dot(A2,R_M.T))    # first equation 6.29
    B1 = np.dot(R_M, A1) - 2 * np.dot(B2,T) # 2nd equation 6.29
    B0 = A0 + np.dot(T.T, (np.dot(B2, T) - \
                           np.dot(R_M, A1))) # 3rd equation 6.29

    return [ B2[0,0], B2[1,1], B2[2,2],
        B2[0,1] + B2[1,0], B2[1,2] + B2[2,1], B2[0,2] + B2[2,0],
        B1[0], B1[1], B1[2], B0]

# factory parameters (input)
p = 10.0
q = 3.0
theta = 0.003

# ellipse parameters
a = (p + q) / 2
b = np.sqrt(p * q) * np.sin(theta)
c = np.sqrt(a**2 - b**2)

# mirror center
yc = (p ** 2 - q ** 2) / 4 / c
zc = -b * np.sqrt(1 - yc ** 2 / a ** 2)
Xc = np.array([0, yc, zc])

# normal to the mirror at center
N = np.array((0, -2 * yc / a ** 2, -2 * zc / b ** 2))
N /= np.sqrt((N**2).sum())

# angle between N and Z
Theta = np.arcsin(N[1])

# rotation matrix
R_M = np.array([[1,0,0],
                [0,np.cos(Theta),-np.sin(Theta)],
                [0,np.sin(Theta),np.cos(Theta)]])

# translation vector
T = -np.dot(R_M,Xc)

# coefficients of the ellipsoid at the centered system
c_in = [1/b**2,1/a**2,1/b**2,0,0,0,0,0,0,-1]

# transformed coeffcients
c_out = rotate_and_translate_coefficients(c_in,R_M,T)

print("coeffs in centered frame: ", c_in)
print("coeffs in local frame: ", c_out)

# results of run
# coeffs in centered frame:
# [3703.715,0.02367,3703.715,0,0,0,0,0,0,-1]
# coeffs in local frame:
# [3703.715,0.0333353,3703.705,0,11.966,0,0,0,-102.564,0]
