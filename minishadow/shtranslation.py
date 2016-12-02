
def shtranslation(a0,qdist1):

# FUNCTION shtranslation,a0,qdist1,file=file
#
# ; qdist: translation vector
#
# IF N_ELEMENTS(qdist1) EQ 0 THEN qdist = [0D,0,0] ELSE qdist=qdist1
#
#
# a0=readsh(a0)
# a=a0

    a = a0.copy()

    a[0,:] = a0[0,:] + qdist1[0]
    a[1,:] = a0[1,:] + qdist1[1]
    a[2,:] = a0[2,:] + qdist1[2]

    return a




