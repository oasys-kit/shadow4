import numpy
from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

if __name__ == "__main__":

    from srxraylib.plot.gol import plot,plot_scatter, set_qt
    set_qt()

    do_plot=True

    # rectangle
    gs = SourceGeometrical(nrays=1000)
    gs.set_spatial_type_rectangle(4,2)
    rays = gs.calculate_rays()

    for i in range(6):
        print("Limits for column %d : %g,%g"%(1+i,rays[:,i].min(),rays[:,i].max()))
    if do_plot:
        plot_scatter(rays[:,0],rays[:,2],xrange=[-3,3],yrange=[-3,3],title="Rectangle")


    # Gaussian
    gs = SourceGeometrical(nrays=1000)
    gs.set_spatial_type_gaussian(2e-6,1e-6)
    rays = gs.calculate_rays()

    for i in range(6):
        print("Std for column %d : %g"%(i+1,rays[:,i].std()))
    if do_plot:
        plot_scatter(1e6*rays[:,0],1e6*rays[:,2],xrange=[-6,6],yrange=[-3,3],title="Gaussian")

    # Ellipse
    gs = SourceGeometrical(nrays=1000)
    gs.set_spatial_type_ellipse(2,1)
    rays = gs.calculate_rays()
    if do_plot:
        plot_scatter(rays[:,0],rays[:,2],xrange=[-3,3],yrange=[-3,3],title="Ellipse")

    # Uniform divergence
    gs = SourceGeometrical(nrays=1000)
    print(">>>>",rays[:,3].shape)
    gs.set_angular_distribution_uniform(-10e-6,5e-6,-3e-6,6e-6)
    rays = gs.calculate_rays()
    if do_plot:
        plot_scatter(1e6*rays[:,3],1e6*rays[:,5],xrange=[-10,10],yrange=[-10,10],title="Uniform divergence")

    # Gaussian divergence
    gs = SourceGeometrical(nrays=5000)
    print(">>>>",rays[:,3].shape)
    gs.set_angular_distribution_gaussian(2e-5,1e-5)

    rays = gs.calculate_rays()
    if do_plot:
        plot_scatter(1e6*rays[:,3],1e6*rays[:,5],title="Gaussian div")

    # Conical divergence
    gs = SourceGeometrical(nrays=5000)
    print(">>>>",rays[:,3].shape)
    gs.set_angular_distribution_cone(2e-5, 1e-5)

    rays = gs.calculate_rays()
    if do_plot:
        plot_scatter(1e6*rays[:,3],1e6*rays[:,5],title="Conical div")


    # energy monochromatic
    gs = SourceGeometrical(nrays=5000)
    # gs.set_energy_distribution_singleline(2000.0,unit='eV')
    # gs.set_energy_distribution_severallines([1000.0,2000.0,3000.,8000],unit='eV')
    # gs.set_energy_distribution_relativeintensities([1000.0,2000.0,3000.,8000],[1.,2.,3,4],unit='eV')
    gs.set_energy_distribution_gaussian(10000.0,2000.0,unit='eV')

    rays = gs.calculate_rays()
    ev = gs._wavenumber_to_energy(rays[:,10])
    if do_plot:
        plot_scatter(rays[:,11],ev,title="Energy Gaussian 10000.0,2000.0")


    # energy external spectrum
    gs = SourceGeometrical(nrays=5000)
    x = numpy.linspace(1000.0,100000,2000)
    y = numpy.exp(- (x-50000)**2 / 2 / 10000**2 )
    # plot(x,y)
    gs.set_energy_distribution_userdefined(x,y,unit='eV')
    rays = gs.calculate_rays()
    ev = gs._wavenumber_to_energy(rays[:,10])
    if do_plot:
        plot_scatter(rays[:,11],ev,title="Energy: sampled from numerical spectrum")

    print(dir(gs))

