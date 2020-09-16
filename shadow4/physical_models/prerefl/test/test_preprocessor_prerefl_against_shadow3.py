import os
from shadow4.physical_models.prerefl.prerefl import PreRefl
import numpy

if __name__ == "__main__":

    # create preprocessor file
    prerefl_file = "Be5_30.dat"
    PreRefl.prerefl(interactive=False, SYMBOL="Be",DENSITY=1.848,FILE=prerefl_file,E_MIN=5000.0,E_MAX=30000.0,E_STEP=100.0)


    #
    #
    # shadow3 fortran
    #
    from shadow4.compatibility.global_definitions import SHADOW3_BINARY
    f = open("prerefl_test.inp",'w')
    f.write("prerefl_test\nBe5_30.dat\n10000")
    f.close()
    print("File written to disk: prerefl_test.inp")
    os.system("%s < prerefl_test.inp > tmp.txt" % SHADOW3_BINARY)
    os.system("cat tmp.txt")


    #
    # python implementation
    #



    a = PreRefl()
    a.read_preprocessor_file(prerefl_file)

    refraction_index = a.get_refraction_index(10000.0,verbose=True)
    print("n = ", refraction_index)

    #
    # compare
    #
    f = open('tmp.txt')  # Open file on read mode
    lines = f.read().split("\n")  # Create a list containing all lines
    f.close()  # Close file


    beta_shadow3 = float(lines[-6][10:])
    delta_shadow3 =  float(lines[-7][10:])
    beta_shadow4 = refraction_index.imag
    delta_shadow4 =  1.0 - refraction_index.real
    print(">>> beta shadow3, shadow4:  ", beta_shadow3, beta_shadow4)
    print(">>> delta shadow3, shadow4: ", delta_shadow3, delta_shadow4 )

    assert (numpy.abs(delta_shadow3 - delta_shadow4) < 1e-8)
    assert (numpy.abs(beta_shadow3 - beta_shadow4) < 1e-11)