###BOOKLISTSTART1###
from amuse.lab import Particles, units 
import pandas as pd
import numpy as np

def new_system():
    '''
    Define the bodies here:
    '''
    particles = Particles(3)
    planet_1 = particles[0]
    planet_1.mass = 1.0 | units.MSun
    planet_1.radius = 1.0 | units.RSun
    planet_1.position = (855251, -804836, -3186) |units.km
    planet_1.velocity = (7.893, 11.894, 0.20642) | (units.m/units.s)

    planet_2 = particles[1]
    planet_2.mass = 0.0025642 | units.MJupiter
    planet_2.radius = 3026.0 | units.km
    planet_2.position = (-0.3767, 0.60159, 0.03930) | units.AU
    planet_2.velocity = (-29.7725, -18.849, 0.795) | units.kms

    planet_3 = particles[2]
    planet_3.mass = 1.0 | units.MEarth
    planet_3.radius = 1.0 | units.REarth
    planet_3.position = (-0.98561, 0.0762, -7.847e-5) | units.AU  
    planet_3.velocity = (-2.927, -29.803, -0.0005327) | units.kms

    particles.move_to_center()
    return particles
###BOOKLISTSTOP1###

###BOOKLISTSTART2###
def integrate_solar_system(particles, end_time):
    '''
    Caclulate the trajectory initial potitions until end_time
    You can use any integrater
    '''
    from amuse.lab import Hermite, nbody_system, Brutus
    convert_nbody = nbody_system.nbody_to_si(particles.mass.sum(),
                                             particles[1].position.length())

    gravity = Hermite(convert_nbody) # you can use any integrator here
    gravity.particles.add_particles(particles)
    planet_1 = gravity.particles[0]
    planet_2 = gravity.particles[1]
    planet_3 = gravity.particles[2]
    def to_float(c):
        return [c[0].value_in(units.AU), c[1].value_in(units.AU), c[2].value_in(units.AU)]
    planet_1_pos = [to_float([planet_1.x,planet_1.y,planet_1.z])]
    planet_2_pos = [to_float([planet_2.x,planet_2.y,planet_2.z])]
    planet_3_pos = [to_float([planet_3.x,planet_3.y,planet_3.z])]
    start_array = np.array([planet_1_pos, planet_2_pos, planet_3_pos]).reshape([-1,3])
    start = pd.DataFrame(data=start_array, columns=["x", "y", "z"])
    
    while gravity.model_time < end_time:
        gravity.evolve_model(gravity.model_time + (1 | units.day)) # define the timesteps you want
        planet_1_pos.append(to_float([planet_1.x,planet_1.y,planet_1.z]))
        planet_2_pos.append(to_float([planet_2.x,planet_2.y,planet_2.z]))
        planet_3_pos.append(to_float([planet_3.x,planet_3.y,planet_3.z]))
    gravity.stop()
    return planet_1_pos, planet_2_pos, planet_3_pos, start
###BOOKLISTSTOP2###

###BOOKLISTSTART3###
def plot_track(planet_1, planet_2, planet_3, start, output_filename):
    '''
    Creates a 3d plot of the orbit
    '''
    import matplotlib.pyplot as plot
    from mpl_toolkits.mplot3d import Axes3D

    x_label = 'x [au]'
    y_label = 'y [au]'
    z_label = 'z [au]'

    fig = plot.figure(figsize=(16, 12), dpi=80)

    colors = ['r','b','g']
    trajectories=[planet_1, planet_2, planet_3]
    ax = fig.add_subplot(1,1,1, projection='3d')
    max_range = 0
    for body, color in zip(trajectories, colors):
        body = np.array(body).reshape([-1,3])
        body = pd.DataFrame(data=body, columns=["x", "y", "z"])
        max_dim = max(max(body["x"]),max(body["y"]),max(body["z"]))
        if max_dim > max_range:
            max_range = max_dim
        ax.plot(body["x"], body["y"], body["z"], color)

    ax.scatter(start["x"], start["y"], start["z"], marker='o', s=200)

    ax.set_xlim([-max_range,max_range])    
    ax.set_ylim([-max_range,max_range])
    ax.set_zlim([-max_range,max_range])

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)

    ax.legend()       

    save_file = output_filename
    plot.savefig(save_file)
    print('\nSaved figure in file', save_file,'\n')
    plot.show()
###BOOKLISTSTOP3###

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-o", 
                      dest="output_filename", default ="SunVenusEarth",
                      help="output filename [%default]")
    return result
    
if __name__ in ('__main__','__plot__'):
    o, arguments  = new_option_parser().parse_args()

    particles = new_system()
    p1, p2, p3, start = integrate_solar_system(particles, 4 | units.yr)
    plot_track(p1, p2, p3, start, o.output_filename)
    